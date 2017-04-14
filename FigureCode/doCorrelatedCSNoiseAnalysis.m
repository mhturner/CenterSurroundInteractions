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

    
    correlationValues = [-1 -0.5 0 0.5 1];
    LinearDeviation = nan(length(populationNodes),length(correlationValues)); %rows are cells, columns correlation values
    OFFcellInds = [];
    ONcellInds = [];
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
        for cc = 1:CorrelatedNoiseNode.children.length
            currentNode = CorrelatedNoiseNode.children(cc);
            currentCorrelationValue = currentNode.splitValue;
            correlationIndex = find(currentCorrelationValue == correlationValues);

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

            %timing stuff:
            timingEpoch = currentNode.childBySplitValue('Center').epochList.firstValue;
            frameRate = timingEpoch.protocolSettings('background:Microdisplay Stage@localhost:monitorRefreshRate');
            if isempty(frameRate)
                frameRate = timingEpoch.protocolSettings('background:Microdisplay_Stage@localhost:monitorRefreshRate');
            end
            FMdata = (riekesuite.getResponseVector(timingEpoch,'Frame Monitor'))';
            [measuredFrameTimes, ~] = getFrameTiming(FMdata,0); %frame flips in data points
            sampleRate = currentNode.epochList.firstValue.protocolSettings('sampleRate'); %Hz
            preTime = currentNode.epochList.firstValue.protocolSettings('preTime') / 1e3; %sec
            stimTime = currentNode.epochList.firstValue.protocolSettings('stimTime') / 1e3; %sec
            preFrames = preTime * frameRate; %frames
            frameDwell = currentNode.epochList.firstValue.protocolSettings('frameDwell'); %frames
            prePoints = preTime * sampleRate; %data points
            stimPoints = stimTime * sampleRate; %data points

            saccadeTimes = measuredFrameTimes(preFrames+1 : frameDwell : end)';
            tempInd = find(saccadeTimes > (prePoints + stimPoints));
            saccadeTimes(tempInd+1 : end) = [];
            
            fixationResponses = nan(3,length(saccadeTimes) - 1);
            for ff = 1:(length(saccadeTimes) - 1)
                tempStart = saccadeTimes(ff);
                tempEnd = saccadeTimes(ff+1);
                fixationResponses(1,ff) = chargeMult * trapz(centerResp.mean(tempStart:tempEnd)) / sampleRate; %pC
                fixationResponses(2,ff) = chargeMult * trapz(surroundResp.mean(tempStart:tempEnd)) / sampleRate; %pC
                fixationResponses(3,ff) = chargeMult * trapz(centerSurroundResp.mean(tempStart:tempEnd)) / sampleRate; %pC
            end
            
            linearDeviation = (fixationResponses(1,:) + fixationResponses(2,:)) - fixationResponses(3,:);
            LinearDeviation(pp,correlationIndex) =  mean(linearDeviation); %#ok<FNDSB>

            if and(CorrelatedNoiseNode.custom.get('isExample'), ismember(currentCorrelationValue, [-1, 0, 1]))
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
        
    end %for cells
    %population plots
    %On cells:
    onMat = LinearDeviation(ONcellInds,:);
    nOn = sum(~isnan(onMat),1);
    meanOn = nanmean(onMat,1);
    errOn = nanstd(onMat,[],1) ./ sqrt(nOn);
    addLineToAxis(correlationValues,meanOn,'meanOn',fig11,[0.7 0.7 0.7],'-','o')
    addLineToAxis(correlationValues,meanOn + errOn,'errOnUp',fig11,[0.7 0.7 0.7],'--','none')
    addLineToAxis(correlationValues,meanOn - errOn,'errOnDown',fig11,[0.7 0.7 0.7],'--','none')
    
    %Off cells:
    offMat = LinearDeviation(OFFcellInds,:);
    nOff = sum(~isnan(offMat),1);
    meanOff = nanmean(offMat,1);
    errOff = nanstd(offMat,[],1) ./ sqrt(nOff);
    addLineToAxis(correlationValues,meanOff,'meanOff',fig11,'k','-','o')
    addLineToAxis(correlationValues,meanOff + errOff,'errOffUp',fig11,'k','--','none')
    addLineToAxis(correlationValues,meanOff - errOff,'errOffDown',fig11,'k','--','none')
    
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