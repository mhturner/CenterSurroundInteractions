function doCorrelatedCSNoiseAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'convertToConductance',true,@islogical);
    
    figDir = '~/Documents/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    convertToConductance = ip.Results.convertToConductance;
    
    figColors = pmkmp(8);

    figure; clf; fig1=gca; %eg Linear filters
    set(fig1,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig1,'XLabel'),'String','Time (s)')
    set(get(fig1,'YLabel'),'String','')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig2=gca; %eg Mean trace: CS and lin sum. Corr = -1
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Time (s)')
    set(get(fig2,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig3=gca; %eg Mean trace: CS and lin sum. Corr = 0
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Time (s)')
    set(get(fig3,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig4=gca; %eg Mean trace: CS and lin sum. Corr = +1
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','Time (s)')
    set(get(fig4,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig5=gca; %eg stim trace: Corr = -1
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','Time (s)')
    set(get(fig5,'YLabel'),'String','Intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig6=gca; %eg stim trace: Corr = 0
    set(fig6,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig6,'XLabel'),'String','Time (s)')
    set(get(fig6,'YLabel'),'String','Intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig7=gca; %eg stim trace: Corr = +1
    set(fig7,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig7,'XLabel'),'String','Time (s)')
    set(get(fig7,'YLabel'),'String','Intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig8=gca; %eg stim cloud: Corr = -1
    set(fig8,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig8,'XLabel'),'String','Center')
    set(get(fig8,'YLabel'),'String','Surround')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig9=gca; %eg stim cloud: Corr = 0
    set(fig9,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig9,'XLabel'),'String','Center')
    set(get(fig9,'YLabel'),'String','Surround')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig10=gca; %eg stim cloud: Corr = +1
    set(fig10,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig10,'XLabel'),'String','Center')
    set(get(fig10,'YLabel'),'String','Surround')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig11=gca; %eg stim cloud: Corr = +1
    set(fig11,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig11,'XLabel'),'String','CS correlation')
    set(get(fig11,'YLabel'),'String','Deviation from linearity')
    set(gcf, 'WindowStyle', 'docked')
    
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
        filterNoiseNode = populationNodes{pp}.childBySplitValue('CenterSurroundNoise').children(1);
        
% % % % % % % % GET LINEAR FILTERS % % % % % % % % % % % % % % % % 
        % center:
        center = getLinearFilterAndPrediction(filterNoiseNode.childBySplitValue('Center').epochList,recType,...
            'seedName','centerNoiseSeed');
        % surround:
        surround = getLinearFilterAndPrediction(filterNoiseNode.childBySplitValue('Surround').epochList,recType,...
            'seedName','surroundNoiseSeed');
        
        if CorrelatedNoiseNode.custom.get('isExample')
            %filters:
            addLineToAxis([center.filterTimeVector],[center.LinearFilter],...
                'center',fig1,figColors(1,:),'-','none')
            addLineToAxis([surround.filterTimeVector],[surround.LinearFilter],...
                'surround',fig1,figColors(4,:),'-','none')
            addLineToAxis(0,0,cellInfo.cellID,fig1,'k','none','none')
        end
        
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
            
            meanLinearDeviation(correlationIndex) = mean(linSum - measuredResponse);
            
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
        figID = ['CScorr_filters_',recID];
        makeAxisStruct(fig1,figID ,'RFSurroundFigs')

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


