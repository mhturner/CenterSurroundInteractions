function doCorrelatedCSNoiseAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'convertToConductance',false,@islogical);
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    convertToConductance = ip.Results.convertToConductance;
    
    figColors = pmkmp(8);
    
    figure; clf; fig2=gca; initFig(fig2,'Time (s)','Response') %eg Mean trace: CS and lin sum. Corr = -1
    figure; clf; fig3=gca; initFig(fig3,'Time (s)','Response') %eg Mean trace: CS and lin sum. Corr = 0
    figure; clf; fig4=gca; initFig(fig4,'Time (s)','Response') %eg Mean trace: CS and lin sum. Corr = +1
    
    figure; clf; fig5=gca; initFig(fig5,'Time (s)','Intensity') %eg stim trace: Corr = -1
    figure; clf; fig6=gca; initFig(fig6,'Time (s)','Intensity') %eg stim trace: Corr = 0
    figure; clf; fig7=gca; initFig(fig7,'Time (s)','Intensity') %eg stim trace: Corr = +1
    
    figure; clf; fig8=gca; initFig(fig8,'Center','Surround') %eg stim cloud: Corr = -1
    figure; clf; fig9=gca; initFig(fig9,'Center','Surround') %eg stim cloud: Corr = 0
    figure; clf; fig10=gca; initFig(fig10,'Center','Surround') %eg stim cloud: Corr = +1
    
    figure; clf; fig11=gca; initFig(fig11,'CS correlation','Deviation from linearity (pC)') %Pop data: dev. from linear vs CS correlation
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(char(node.descendentsDepthFirst(nn).splitKey),...
                '@(list)splitOnShortProtocolID(list)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    meanDiff.shuffle = [];
    meanDiff.control = [];
    ONcellInds = [];
    OFFcellInds = [];
    for pp = 1:length(populationNodes)
        cellInfo = getCellInfoFromEpochList(populationNodes{pp}.epochList);
        recType = getRecordingTypeFromEpochList(populationNodes{pp}.epochList);
        if (convertToConductance)
            recType = [recType, ', conductance']; %#ok<AGROW>
        end
        if strcmp(cellInfo.cellType,'ONparasol')
            ONcellInds = cat(2,ONcellInds,pp);
            plotColor = [0.7 0.7 0.7];
        elseif strcmp(cellInfo.cellType,'OFFparasol')
            OFFcellInds = cat(2,OFFcellInds,pp);
            plotColor = 'k';
        end
        
        CorrelatedNoiseNode = populationNodes{pp}.childBySplitValue('CorrelatedCSNoise');
        
% % % % % % % % DO ADDITIVITY ANALYSIS % % % % % % % % % % % % % % % %
        meanLinearDeviation = [];
        correlationValues = [];
        for correlationIndex = 1:CorrelatedNoiseNode.children.length
            currentNode = CorrelatedNoiseNode.children(correlationIndex);
            currentCorrelationValue = currentNode.splitValue;
            correlationValues(correlationIndex) = currentCorrelationValue;

            centerResp = getMeanResponseTrace(currentNode.childBySplitValue('Center').epochList,recType);
            surroundResp = getMeanResponseTrace(currentNode.childBySplitValue('Surround').epochList,recType);
            centerSurroundResp = getMeanResponseTrace(currentNode.childBySplitValue('Center-Surround').epochList,recType);
            
            timeVec = centerSurroundResp.timeVector;
            measuredResponse = centerSurroundResp.mean;
            linSum = centerResp.mean + surroundResp.mean;
            
            if strcmp(recType,'exc')
                chargeMult = -1;
            else
                chargeMult = 1;
            end

            meanLinearDeviation(correlationIndex) = chargeMult * mean(linSum - measuredResponse);
            
            %reconstruct noise stimuli:
            backgroundIntensity = currentNode.epochList.firstValue.protocolSettings('backgroundIntensity');
            sampleRate = currentNode.epochList.firstValue.protocolSettings('sampleRate');
            prePts = sampleRate * currentNode.epochList.firstValue.protocolSettings('preTime') / 1e3;
            frameDwell = currentNode.epochList.firstValue.protocolSettings('frameDwell');
            lightCrafterFlag = currentNode.epochList.firstValue.protocolSettings.keySet.contains('background:LightCrafter Stage@localhost:lightCrafterPatternRate');
            if ~lightCrafterFlag
                lightCrafterFlag = currentNode.epochList.firstValue.protocolSettings.keySet.contains('background:LightCrafter_Stage@localhost:lightCrafterPatternRate');
            end
            FMdata = (riekesuite.getResponseVector(currentNode.epochList.firstValue,'Frame Monitor'))';
            frameTimes = getFrameTiming(FMdata,lightCrafterFlag);
            updateLength = frameDwell*mean(diff(frameTimes));
            
            centerStim = convertJavaArrayList(currentNode.epochList.firstValue.protocolSettings('centerNoiseArray'));
            surroundStim = convertJavaArrayList(currentNode.epochList.firstValue.protocolSettings('surroundNoiseArray'));
            
            centerStim_pt = centerStim;
            surroundStim_pt = surroundStim;
            
            centerStim = [backgroundIntensity.* ones(1,prePts), kron(centerStim,ones(1,round(updateLength))), backgroundIntensity.* ones(1,prePts)];
            surroundStim = [backgroundIntensity.* ones(1,prePts), kron(surroundStim,ones(1,round(updateLength))), backgroundIntensity.* ones(1,prePts)];
            stimTimeVec = (1:length(centerStim)) / sampleRate;

            if and(CorrelatedNoiseNode.custom.get('isExample'), ismember(currentCorrelationValue, [-1, 0, 1]))
                if currentCorrelationValue == -1
                    figA = fig2; figB = fig5; figC = fig8;
                elseif currentCorrelationValue == 0
                    figA = fig3; figB = fig6; figC = fig9;
                elseif currentCorrelationValue == 1
                    figA = fig4; figB = fig7; figC = fig10;
                end
                addLineToAxis(timeVec,measuredResponse,'measured',figA,'k','-','none')
                addLineToAxis(timeVec,linSum,'linSum',figA,[0.7 0.7 0.7],'-','none')
                
                addLineToAxis(stimTimeVec,centerStim,'center',figB,figColors(1,:),'-','none')
                addLineToAxis(stimTimeVec,surroundStim,'surround',figB,figColors(4,:),'-','none')
                
                addLineToAxis(centerStim_pt,surroundStim_pt,'data',figC,'k','none','o')
            end %eg cell plots
        end %for correlation values
        
        addLineToAxis(correlationValues,meanLinearDeviation,['data', num2str(pp)],fig11,plotColor,'-','o')

    end %for cells

    recID = getRecordingTypeFromEpochList(currentNode.epochList);
    if (exportFigs)
        figID = ['CScorr_respNeg_',recID];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')

        figID = ['CScorr_respZero_',recID];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')

        figID = ['CScorr_respPos_',recID];
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')

        figID = ['CScorr_stimNeg_',recID];
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')
        
        figID = ['CScorr_stimZero_',recID];
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')
        
        figID = ['CScorr_stimPos_',recID];
        makeAxisStruct(fig7,figID ,'RFSurroundFigs')
        
        figID = ['CScorr_cloudNeg_',recID];
        makeAxisStruct(fig8,figID ,'RFSurroundFigs')
        
        figID = ['CScorr_cloudZero_',recID];
        makeAxisStruct(fig9,figID ,'RFSurroundFigs')
        
        figID = ['CScorr_cloudPos_',recID];
        makeAxisStruct(fig10,figID ,'RFSurroundFigs')
        
        figID = ['CScorr_popDeviation_',recID];
        makeAxisStruct(fig11,figID ,'RFSurroundFigs')
        
    end
end