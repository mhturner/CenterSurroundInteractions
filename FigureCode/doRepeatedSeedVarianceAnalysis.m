function doRepeatedSeedVarianceAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));

    ip.parse(node,varargin{:});
    node = ip.Results.node;

    figure; clf; fig1=gca; initFig(fig1,'Time (s)','') % Linear filters
    figure; clf; fig2=gca; initFig(fig2,'Linear prediction (nS)','Measured (nS)') % Nonlinearities
    figure; clf; fig3=gca; initFig(fig3,'Linear prediction (nS)','Measured (nS)') % Independent nonlinearities
    figure; clf; fig4=gca; initFig(fig4,'Linear prediction (nS)','Measured (nS)') % Shared nonlinearity

    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(currentStimulus)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
        
    ONcellInds = [];
    OFFcellInds = [];
    reliabilityVals = [];
    for pp = 1:length(populationNodes)
        currentNode = populationNodes{pp};

        cellInfo = getCellInfoFromEpochList(currentNode.epochList);

        if strcmp(cellInfo.cellType,'ONparasol')
            ONcellInds = cat(2,ONcellInds,pp);
        elseif strcmp(cellInfo.cellType,'OFFparasol')
            OFFcellInds = cat(2,OFFcellInds,pp);
        end
        CSnode = currentNode.childBySplitValue('Center-Surround');

        centerSurround = getLinearFilterAndPrediction(CSnode.epochList,'exc',...
            'seedName','centerNoiseSeed','numberOfBins',10);
        
        epochLen = length(centerSurround.measuredResponse) / centerSurround.n;
        meanTrace = mean(reshape(centerSurround.measuredResponse,centerSurround.n,epochLen),1);
        
        r2vals = [];
        for testingEpoch = 1:(centerSurround.n)
            testDataInds = ((testingEpoch-1)*epochLen + 1):(testingEpoch*epochLen);
            ss_total = sum((centerSurround.measuredResponse-mean(centerSurround.measuredResponse)).^2);
            ss_resid = sum((meanTrace-centerSurround.measuredResponse(testDataInds)).^2);
            r2Vals(testingEpoch) = 1-ss_resid/ss_total;
        end
        
        reliabilityVals(pp) = mean(r2Vals);
        
    end
  
    disp(mean(reliabilityVals(ONcellInds)))
    disp(mean(reliabilityVals(OFFcellInds)))
end