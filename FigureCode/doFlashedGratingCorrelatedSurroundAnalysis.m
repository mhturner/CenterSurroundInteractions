function doFlashedGratingCorrelatedSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    addParameter(ip,'exportFigs',true,@islogical);

    ip.parse(node,varargin{:});
    node = ip.Results.node;
    metric = ip.Results.metric;
    exportFigs = ip.Results.exportFigs;
    
    
    targetContrast = 0.5;
    targetIntensityValues = [0.3 0.4 0.5 0.6 0.7];
    
    figure; clf; fig1=gca; initFig(fig1,'Center Intensity','NLI') % NLI vs. int for diff surrounds
    
    figure; clf; fig2=gca; initFig(fig2,'Center Intensity','NLI') % MEAN NLI vs. int for diff surrounds
 
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(gratingContrast)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    nliMat = [];
    for pp = 1:length(populationNodes)
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        
        cellNode.children.length
        gratingNode = cellNode.childBySplitValue(0.5);
            
        respMat = nan(length(targetIntensityValues), 3, 2); %intensity x surroundType x image/disc;
        for ii = 1:length(targetIntensityValues)%for each currentIntensity
            
            intensityNode = gratingNode.childBySplitValue(targetIntensityValues(ii));

            for ss = 1:3 %for each surround type (acorr, corr, none)
                surroundNode = intensityNode.children(ss);

                imageResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('image').epochList,recType);
                discResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('intensity').epochList,recType);

                respMat(ii,ss,1) = imageResp.(metric).mean;
                respMat(ii,ss,2) = discResp.(metric).mean;

            end %end for surroundTag
        end %end for currentIntensity      
        
        nliMat(:,:,pp) = (respMat(:,:,1) - respMat(:,:,2)) ./ (respMat(:,:,1) + respMat(:,:,2));
        
        addLineToAxis(targetIntensityValues,nliMat(:,1,pp),'acorr',fig1,'r','-','o')
        addLineToAxis(targetIntensityValues,nliMat(:,2,pp),'corr',fig1,'g','-','o')
        addLineToAxis(targetIntensityValues,nliMat(:,3,pp),'none',fig1,'k','-','o')
    end %end for cell

    
    addLineToAxis(targetIntensityValues,nanmean(nliMat(:,1,:),3),'acorr',fig2,'r','-','o')
    addLineToAxis(targetIntensityValues,nanmean(nliMat(:,2,:),3),'corr',fig2,'g','-','o')
    addLineToAxis(targetIntensityValues,nanmean(nliMat(:,3,:),3),'none',fig2,'k','-','o')
    if (exportFigs)

    end
end