function doCorrelatedCSadditivityAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'convertToConductance',false,@islogical);
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    convertToConductance = ip.Results.convertToConductance;
    
    figColors = pmkmp(8);
    
    figure; clf; fig2=gca; initFig(fig2,'C/S','C+S') %eg Mean trace: CS and lin sum. Corr = -1
    figure(30); clf;
   
    targetCorrelation = 0;
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(char(node.descendentsDepthFirst(nn).splitKey),...
                '@(list)splitOnShortProtocolID(list)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    
    for pp = 1:length(populationNodes)
        cellInfo = getCellInfoFromEpochList(populationNodes{pp}.epochList);
        recType = getRecordingTypeFromEpochList(populationNodes{pp}.epochList);
        if (convertToConductance)
            recType = [recType, ', conductance']; %#ok<AGROW>
        end
        CorrelatedNoiseNode = populationNodes{pp}.childBySplitValue('CorrelatedCSNoise');
        currentNode = CorrelatedNoiseNode.childBySplitValue(targetCorrelation);
% % % % % % % % DO ADDITIVITY ANALYSIS % % % % % % % % % % % % % % % %
        centerResp = getMeanResponseTrace(currentNode.childBySplitValue('Center').epochList,recType);
        surroundResp = getMeanResponseTrace(currentNode.childBySplitValue('Surround').epochList,recType);
        centerSurroundResp = getMeanResponseTrace(currentNode.childBySplitValue('Center-Surround').epochList,recType);

        timeVec = centerSurroundResp.timeVector;
        measuredResponse = centerSurroundResp.mean;
        linSum = centerResp.mean + surroundResp.mean;
        
        figure(30);
        fh = subplot(3,3,pp);
        binAndPlotPopulationData(measuredResponse,linSum,50,fh,'b')
        limDown = min([measuredResponse,linSum]); limUp = max([measuredResponse,linSum]);
        addLineToAxis([limDown,limUp],[limDown,limUp],'unity',fh,'k','--','none');
        if pp == 2
            title(['Corr = ',num2str(targetCorrelation)])
        elseif pp == length(populationNodes)
            xlabel('C/S'); ylabel('C + S')
        end
        set(fh,'FontSize',10)
        
    end %for cells
    
    recID = getRecordingTypeFromEpochList(currentNode.epochList);
    if (exportFigs)
% %         figID = ['CSadd_eg_',recID];
% %         makeAxisStruct(fig2,figID ,'RFSurroundFigs')

        
    end
end