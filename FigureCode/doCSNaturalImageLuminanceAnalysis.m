function doCSNaturalImageLuminanceAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    
    figDir = '~/Dropbox/RiekeLab/Analysis/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    
    figColors = pmkmp(8);
    
    figure; clf; fig2=gca; initFig(fig2,'Time (s)','Response') %eg Mean trace: C, S, indep, no shuffle
    figure; clf; fig3=gca; initFig(fig3,'Time (s)','Response') %eg Mean trace: CS together and linear sum, no shuffle
    
    figure; clf; fig4=gca; initFig(fig4,'Time (s)','Response') %eg Mean trace: C, S, indep, with shuffle
    figure; clf; fig5=gca; initFig(fig5,'Time (s)','Response') %eg Mean trace: CS together and linear sum, with shuffle
    
    figure; clf; fig6=gca; initFig(fig6,'R(C+S)','R(C) + R(S)') %eg scatter plot of R(C) + R(S) vs R(C+S)
    figure; clf; fig7=gca; initFig(fig7,'Control Diff. (pC)','Shuffle Diff. (pC)') %population. Difference (linear sum - measured) for control vs shuffled. ON and OFF cells

    figure; clf; fig8=gca; initFig(fig8,'Time(s)','Intensity') %eg. C and S stims for control
    figure; clf; fig9=gca; initFig(fig9,'Time(s)','Intensity') %eg. C and S stims for shuffled
    
    figure; clf; fig10=gca; initFig(fig10,'Center intensity','Surround intensity') %eg. C and S stim scatter, control
    figure; clf; fig11=gca; initFig(fig11,'Center intensity','Surround intensity') %eg. C and S stim scatter, shuffled
    
    figure; clf; fig12=gca; initFig(fig12,'Intensity','Probability') %eg. Nat image intensity histogram
    
    figure; clf; fig13=gca; initFig(fig13,'Sparsity natural CS','Sparsity shuffled CS') %sparsity CS corr vs shuffled
    
    figure(30); clf; 
    
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
    sparsity = struct;
    ONcellInds = [];
    OFFcellInds = [];
    for pp = 1:length(populationNodes)
        cellInfo = getCellInfoFromEpochList(populationNodes{pp}.epochList);
        recType = getRecordingTypeFromEpochList(populationNodes{pp}.epochList);
        if strcmp(cellInfo.cellType,'ONparasol')
            ONcellInds = cat(2,ONcellInds,pp);
        elseif strcmp(cellInfo.cellType,'OFFparasol')
            OFFcellInds = cat(2,OFFcellInds,pp);
        end
        
        ImageLuminanceNode = populationNodes{pp}.childBySplitValue('CSNaturalImageLuminance');

% % % % % % % % DO ADDITIVITY ANALYSIS % % % % % % % % % % % % % % % %
        fixationResponses = nan(6,1); %rows are controlCenter, controlSurround,...
        fixationCount = 0;
        frameCount = 0;
        for imageIndex = 1:ImageLuminanceNode.children.length
            currentNode = ImageLuminanceNode.children(imageIndex);
            for ss = 1:2
                if currentNode.children(ss).splitValue == 0 %no shuffle CS
                    controlCSNode = currentNode.children(ss);
                elseif currentNode.children(ss).splitValue == 1 %shuffle CS
                    shuffleCSNode = currentNode.children(ss);
                end
            end
            
            if strcmp(recType,'extracellular')
                attachSpikeBinary = true;
            else
                attachSpikeBinary = false;
            end

            controlCenter = getMeanResponseTrace(controlCSNode.childBySplitValue('Center').epochList,recType,'attachSpikeBinary',attachSpikeBinary);
            controlSurround = getMeanResponseTrace(controlCSNode.childBySplitValue('Surround').epochList,recType,'attachSpikeBinary',attachSpikeBinary);
            controlCenterSurround = getMeanResponseTrace(controlCSNode.childBySplitValue('Center-Surround').epochList,recType,'attachSpikeBinary',attachSpikeBinary);

            shuffleCenter = getMeanResponseTrace(shuffleCSNode.childBySplitValue('Center').epochList,recType,'attachSpikeBinary',attachSpikeBinary);
            shuffleSurround = getMeanResponseTrace(shuffleCSNode.childBySplitValue('Surround').epochList,recType,'attachSpikeBinary',attachSpikeBinary);
            shuffleCenterSurround = getMeanResponseTrace(shuffleCSNode.childBySplitValue('Center-Surround').epochList,recType,'attachSpikeBinary',attachSpikeBinary);

            figure(30);
            fh = subplot(3,3,pp);
            measuredResponse = controlCenterSurround.mean;
            linSum = controlCenter.mean + controlSurround.mean;
            binAndPlotPopulationData(measuredResponse,linSum,50,fh,'b')
            limDown = min([measuredResponse,linSum]); limUp = max([measuredResponse,linSum]);
            addLineToAxis([limDown,limUp],[limDown,limUp],'unity',fh,'k','--','none');
            if pp == length(populationNodes)
                xlabel('C/S'); ylabel('C + S')
            end
            set(fh,'FontSize',10)
            
            
            %timing stuff:
            timingEpoch = controlCSNode.childBySplitValue('Center').epochList.firstValue;
            frameRate = timingEpoch.protocolSettings('background:Microdisplay Stage@localhost:monitorRefreshRate');
            if isempty(frameRate)
                frameRate = timingEpoch.protocolSettings('background:Microdisplay_Stage@localhost:monitorRefreshRate');
            end
            FMdata = (riekesuite.getResponseVector(timingEpoch,'Frame Monitor'))';
            [measuredFrameTimes, ~] = getFrameTiming(FMdata,0);
            
            sampleRate = currentNode.epochList.firstValue.protocolSettings('sampleRate');
            preTime = currentNode.epochList.firstValue.protocolSettings('preTime') / 1e3; %sec
            stimTime = currentNode.epochList.firstValue.protocolSettings('stimTime') / 1e3; %sec
            tailTime = currentNode.epochList.firstValue.protocolSettings('tailTime') / 1e3; %sec
            totalTime = preTime + stimTime + tailTime;
            fixationDuration = currentNode.epochList.firstValue.protocolSettings('fixationDuration') / 1e3; %sec
            targetSaccadeTimes = 0:fixationDuration:(stimTime-fixationDuration);
            
            idealFrameTimes = 0:(1/frameRate):totalTime;
            stateTime = idealFrameTimes - preTime;
            saccadeTimes = nan(size(targetSaccadeTimes));
            for ff = 1:length(targetSaccadeTimes)
                tempInd = find(stateTime >= targetSaccadeTimes(ff),1);
                saccadeTimes(ff) = measuredFrameTimes(tempInd);
            end
            saccadeTimes(end+1) = saccadeTimes(end) + mean(diff(saccadeTimes));
            saccadeTimes = round(saccadeTimes);
            noFixations = currentNode.epochList.firstValue.protocolSettings('stimTime') / ...
                currentNode.epochList.firstValue.protocolSettings('fixationDuration');
  
            % get fixation responses for this image
            for ff = 1:noFixations
                fixationCount = fixationCount + 1;
                tempStart = saccadeTimes(ff);
                tempEnd = saccadeTimes(ff+1);
                
                if strcmp(recType,'extracellular') %spikes: spike count
                    %count spikes in each trial on this fixation. Mean
                    %across trials:
                    fixationResponses(1,fixationCount) = mean(sum(controlCenter.binary(:,tempStart:tempEnd),2),1);
                    fixationResponses(2,fixationCount) = mean(sum(controlSurround.binary(:,tempStart:tempEnd),2),1); 
                    fixationResponses(3,fixationCount) = mean(sum(controlCenterSurround.binary(:,tempStart:tempEnd),2),1); 
                    
                    fixationResponses(4,fixationCount) = mean(sum(shuffleCenter.binary(:,tempStart:tempEnd),2),1); 
                    fixationResponses(5,fixationCount) = mean(sum(shuffleSurround.binary(:,tempStart:tempEnd),2),1); 
                    fixationResponses(6,fixationCount) = mean(sum(shuffleCenterSurround.binary(:,tempStart:tempEnd),2),1); 
                    
                else %currents: charge transfer
                    if strcmp(recType,'exc')
                        chargeMult = -1;
                    else
                        chargeMult = 1;
                    end
                    fixationResponses(1,fixationCount) = chargeMult * trapz(controlCenter.mean(tempStart:tempEnd)) / sampleRate; %pC
                    fixationResponses(2,fixationCount) = chargeMult * trapz(controlSurround.mean(tempStart:tempEnd)) / sampleRate;
                    fixationResponses(3,fixationCount) = chargeMult * trapz(controlCenterSurround.mean(tempStart:tempEnd)) / sampleRate;

                    fixationResponses(4,fixationCount) = chargeMult * trapz(shuffleCenter.mean(tempStart:tempEnd)) / sampleRate;
                    fixationResponses(5,fixationCount) = chargeMult * trapz(shuffleSurround.mean(tempStart:tempEnd)) / sampleRate;
                    fixationResponses(6,fixationCount) = chargeMult * trapz(shuffleCenterSurround.mean(tempStart:tempEnd)) / sampleRate;
                end
            end %end for fixations

            if strcmp(recType,'extracellular')
                [frameTimes, ~] = getFrameTiming(FMdata,0);
                frameTimes = frameTimes((preFrames + 1) : (preFrames + stimFrames + 1));
                for ff = 1:(length(frameTimes)-1)
                    frameCount = frameCount + 1;
                    startPt = frameTimes(ff);
                    endPt = frameTimes(ff+1);
                    controlCenter.frameResp(frameCount) = sum(mean(controlCenter.binary(:,startPt:endPt),1)); %spikes
                    controlSurround.frameResp(frameCount) = sum(mean(controlSurround.binary(:,startPt:endPt),1)); %spikes
                    controlCenterSurround.frameResp(frameCount) = sum(mean(controlCenterSurround.binary(:,startPt:endPt),1)); %spikes

                    shuffleCenter.frameResp(frameCount) = sum(mean(shuffleCenter.binary(:,startPt:endPt),1)); %spikes
                    shuffleSurround.frameResp(frameCount) = sum(mean(shuffleSurround.binary(:,startPt:endPt),1)); %spikes
                    shuffleCenterSurround.frameResp(frameCount) = sum(mean(shuffleCenterSurround.binary(:,startPt:endPt),1)); %spikes
                end
            end
            

            if currentNode.custom.get('isExample')
                %control:
                %   indep
                addLineToAxis(controlCenter.timeVector,controlCenter.mean,...
                    'center',fig2,figColors(1,:),'-','none')
                addLineToAxis(controlSurround.timeVector,controlSurround.mean,...
                    'surround',fig2,figColors(4,:),'-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig2,'k','none','none')
                %   lin sum vs combined
                addLineToAxis(controlCenter.timeVector,controlCenter.mean + controlSurround.mean,...
                    'linSum',fig3,[0.7 0.7 0.7],'-','none')
                addLineToAxis(controlCenterSurround.timeVector,controlCenterSurround.mean,...
                    'measuredCS',fig3,'k','-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig3,'k','none','none')
                
                %shuffle
                %   indep
                addLineToAxis(shuffleCenter.timeVector,shuffleCenter.mean,...
                    'center',fig4,figColors(1,:),'-','none')
                addLineToAxis(shuffleSurround.timeVector,shuffleSurround.mean,...
                    'surround',fig4,figColors(4,:),'-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig4,'k','none','none')
                %   lin sum vs combined
                addLineToAxis(shuffleCenter.timeVector,shuffleCenter.mean + shuffleSurround.mean,...
                    'linSum',fig5,[0.7 0.7 0.7],'-','none')
                addLineToAxis(shuffleCenterSurround.timeVector,shuffleCenterSurround.mean,...
                    'measuredCS',fig5,'k','-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig5,'k','none','none')
                
                % scatter plot with response per fixation
                % just for this example image
                egIndsToPull = (size(fixationResponses,2) - noFixations + 1):size(fixationResponses,2);
                addLineToAxis(fixationResponses(3,egIndsToPull),...
                    fixationResponses(1,egIndsToPull)+fixationResponses(2,egIndsToPull),...
                    'control',fig6,'k','none','o')
                addLineToAxis(fixationResponses(6,egIndsToPull),...
                    fixationResponses(4,egIndsToPull)+fixationResponses(5,egIndsToPull),...
                    'shuff',fig6,'g','none','o')
                
                    %e.g. cell stats:
                [~,p] = ttest(fixationResponses(3,egIndsToPull),...
                    fixationResponses(1,egIndsToPull)+fixationResponses(2,egIndsToPull));
                disp('e.g. cell meas. vs. linSum, Control:')
                disp(['p = ',num2str(p)])
                
                [~,p] = ttest(fixationResponses(6,egIndsToPull),...
                    fixationResponses(4,egIndsToPull)+fixationResponses(5,egIndsToPull));
                disp('e.g. cell meas. vs. linSum, Shuffled:')
                disp(['p = ',num2str(p)])
                
                downLim = min([fixationResponses(6,egIndsToPull), fixationResponses(1,egIndsToPull)]);
                upLim = max([fixationResponses(6,egIndsToPull), fixationResponses(1,egIndsToPull)]);
                addLineToAxis([downLim upLim],[downLim upLim],'unity',fig6,'k','--','none')
                meanControl = [mean(fixationResponses(3,egIndsToPull)),...
                    mean(fixationResponses(1,egIndsToPull)+fixationResponses(2,egIndsToPull))];
                semControl = [std(fixationResponses(3,egIndsToPull)) / sqrt(length(egIndsToPull)),...
                    std(fixationResponses(1,egIndsToPull)+fixationResponses(2,egIndsToPull)) / sqrt(length(egIndsToPull))];
                
                meanShuffle = [mean(fixationResponses(6,egIndsToPull)),...
                    mean(fixationResponses(4,egIndsToPull)+fixationResponses(5,egIndsToPull))];
                semShuffle = [std(fixationResponses(6,egIndsToPull)) / sqrt(length(egIndsToPull)),...
                    std(fixationResponses(4,egIndsToPull)+fixationResponses(5,egIndsToPull)) / sqrt(length(egIndsToPull))];
                
                addLineToAxis(meanControl(1),meanControl(2),'meanControl',fig6,'k','none','x')
                addLineToAxis([meanControl(1)-semControl(1), meanControl(1)+semControl(1)],...
                    [meanControl(2), meanControl(2)],'errControlX',fig6,'k','-','none')
                addLineToAxis([meanControl(1), meanControl(1)],...
                    [meanControl(2)-semControl(2), meanControl(2)+semControl(2)],'errControlY',fig6,'k','-','none')
                
                addLineToAxis(meanShuffle(1),meanShuffle(2),'meanShuffle',fig6,'r','none','x')
                addLineToAxis([meanShuffle(1)-semShuffle(1), meanShuffle(1)+semShuffle(1)],...
                    [meanShuffle(2), meanShuffle(2)],'errShuffleX',fig6,'r','-','none')
                addLineToAxis([meanShuffle(1), meanShuffle(1)],...
                    [meanShuffle(2)-semShuffle(2), meanShuffle(2)+semShuffle(2)],'errShuffleY',fig6,'r','-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig6,'k','none','none')
                % c and s stims
                %   traces
                startPoint = saccadeTimes(1);
                fixationPoints = fixationDuration * sampleRate;
                pad = ones(1,startPoint) .* controlCSNode.epochList.firstValue.protocolSettings('backgroundIntensity');
                controlCStim_fix = convertJavaArrayList(controlCSNode.epochList.firstValue.protocolSettings('CenterIntensity'));
                controlSStim_fix = convertJavaArrayList(controlCSNode.epochList.firstValue.protocolSettings('SurroundIntensity'));
                controlCStim = [pad, kron(controlCStim_fix,ones(1,fixationPoints)), pad];
                controlSStim = [pad, kron(controlSStim_fix,ones(1,fixationPoints)), pad];

                shuffleCStim_fix = convertJavaArrayList(shuffleCSNode.epochList.firstValue.protocolSettings('CenterIntensity'));
                shuffleSStim_fix = convertJavaArrayList(shuffleCSNode.epochList.firstValue.protocolSettings('SurroundIntensity'));
                shuffleCStim = [pad, kron(shuffleCStim_fix,ones(1,fixationPoints)), pad];
                shuffleSStim = [pad, kron(shuffleSStim_fix,ones(1,fixationPoints)), pad];
                timeVec = (0:(length(controlCStim)-1)) ./ sampleRate;

                addLineToAxis(timeVec,controlCStim,'C',fig8,figColors(1,:),'-','none')
                addLineToAxis(timeVec,controlSStim,'S',fig8,figColors(4,:),'-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig8,'k','none','none')
                
                addLineToAxis(timeVec,shuffleCStim,'C',fig9,figColors(1,:),'-','none')
                addLineToAxis(timeVec,shuffleSStim,'S',fig9,figColors(4,:),'-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig9,'k','none','none')
                
                %   scatter plot of c,s fixations
                pRes = polyfit(controlCStim_fix,controlSStim_fix,1);
                downFit = min([controlCStim_fix, controlSStim_fix]);
                upFit = max([controlCStim_fix, controlSStim_fix]);
                rho = corr(controlCStim_fix',controlSStim_fix');
                addLineToAxis(controlCStim_fix,controlSStim_fix,'pts',fig10,'k','none','o')
                addLineToAxis([downFit upFit],polyval(pRes,[downFit upFit]),['fit',num2str(rho)],fig10,'k','-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig10,'k','none','none')
                
                pRes = polyfit(shuffleCStim_fix,shuffleSStim_fix,1);
                downFit = min([shuffleCStim_fix, shuffleSStim_fix]);
                upFit = max([shuffleCStim_fix, shuffleSStim_fix]);
                rho = corr(shuffleCStim_fix',shuffleSStim_fix');
                addLineToAxis(shuffleCStim_fix,shuffleSStim_fix,'pts',fig11,'r','none','o')
                addLineToAxis([downFit upFit],polyval(pRes,[downFit upFit]),['fit',num2str(rho)],fig11,'r','-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig11,'k','none','none')
                
                %natural image and intensity histogram
                resourcesDir = '~/Dropbox/RiekeLab/Analysis/MATLAB/turner-package/resources/';
                stimSet = 'VHsubsample_20160105';
                load([resourcesDir, controlCSNode.epochList.firstValue.protocolSettings('currentStimSet')]);
                imageName = luminanceData(controlCSNode.epochList.firstValue.protocolSettings('imageIndex')).ImageName;
                fileId=fopen([resourcesDir, stimSet, '/', imageName,'.iml'],'rb','ieee-be');
                img = fread(fileId, [1536,1024], 'uint16');
                img = double(img);
                img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
                fh = figure(21); clf;
                imagesc(img'); colormap(gray); axis image; axis off;
                drawnow;
                figID = 'egNatImage';
                print(fh,[figDir,figID],'-depsc')
                
                [nn, cc] = histcounts(img(:),50,'Normalization','probability');
                binCtrs = cc(1:end-1) + mean(diff(cc));
                addLineToAxis(binCtrs,nn,'hist',fig12,'k','-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig12,'k','none','none')
                
            end %eg cell plots
        end %for images
        
        [~, sparsity.controlCenter(pp)] = getActivityRatio(fixationResponses(1,:));
        [~, sparsity.controlSurround(pp)] = getActivityRatio(fixationResponses(2,:));
        [~, sparsity.controlCenterSurround(pp)] = getActivityRatio(fixationResponses(3,:));
        [~, sparsity.controlLinSum(pp)] = getActivityRatio(fixationResponses(1,:) + fixationResponses(2,:));

        [~, sparsity.shuffleCenter(pp)] = getActivityRatio(fixationResponses(4,:));
        [~, sparsity.shuffleSurround(pp)] = getActivityRatio(fixationResponses(5,:));
        [~, sparsity.shuffleCenterSurround(pp)] = getActivityRatio(fixationResponses(6,:));
        [~, sparsity.shuffleLinSum(pp)] = getActivityRatio(fixationResponses(4,:) + fixationResponses(5,:));

        meanDiff.control(pp) = mean((fixationResponses(1,:) + fixationResponses(2,:)) - fixationResponses(3,:));
        meanDiff.shuffle(pp) = mean((fixationResponses(4,:) + fixationResponses(5,:)) - fixationResponses(6,:));
    end
    
    %Population data:
    %   On...
    addLineToAxis(meanDiff.control(ONcellInds),meanDiff.shuffle(ONcellInds),...
        'ONcells',fig7,'b','none','o')
    meanON_x = mean(meanDiff.control(ONcellInds));
    meanON_y = mean(meanDiff.shuffle(ONcellInds));
    errON_x = std(meanDiff.control(ONcellInds)) ./ sqrt(length(ONcellInds));
    errON_y = std(meanDiff.shuffle(ONcellInds)) ./ sqrt(length(ONcellInds));
    addLineToAxis(meanON_x,meanON_y,'ONmean',fig7,'b','none','.')
    addLineToAxis(meanON_x + [errON_x, -errON_x],[meanON_y meanON_y],'ONerr_x',fig7,'b','-','none')
    addLineToAxis([meanON_x meanON_x],meanON_y + [errON_y, -errON_y],'ONerr_y',fig7,'b','-','none')
    %       pop stats:
    [~,p] = ttest(meanDiff.control(ONcellInds),meanDiff.shuffle(ONcellInds));
    disp('Control vs. shuffled On:')
    disp(['p = ',num2str(p)])
    disp(['n = ',num2str(length(ONcellInds))])
    
    %   Off...
    addLineToAxis(meanDiff.control(OFFcellInds),meanDiff.shuffle(OFFcellInds),...
        'OFFcells',fig7,'r','none','o')
    meanOFF_x = mean(meanDiff.control(OFFcellInds));
    meanOFF_y = mean(meanDiff.shuffle(OFFcellInds));
    errOFF_x = std(meanDiff.control(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    errOFF_y = std(meanDiff.shuffle(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    addLineToAxis(meanOFF_x,meanOFF_y,'OFFmean',fig7,'r','none','.')
    addLineToAxis(meanOFF_x + [errOFF_x, -errOFF_x],[meanOFF_y meanOFF_y],'OFFerr_x',fig7,'r','-','none')
    addLineToAxis([meanOFF_x meanOFF_x],meanOFF_y + [errOFF_y, -errOFF_y],'OFFerr_y',fig7,'r','-','none')
    %       pop stats:
    [~,p] = ttest(meanDiff.control(OFFcellInds),meanDiff.shuffle(OFFcellInds));
    disp('Control vs. shuffled Off:')
    disp(['p = ',num2str(p)])
    disp(['n = ',num2str(length(OFFcellInds))])
    
    downLim = min([meanDiff.control, meanDiff.shuffle, 0]);
    upLim = 1.1*max([meanDiff.control, meanDiff.shuffle]);
    addLineToAxis([downLim upLim],[downLim upLim],'unity',fig7,'k','--','none')
    
    addLineToAxis(sparsity.controlCenterSurround(ONcellInds),sparsity.shuffleCenterSurround(ONcellInds),...
        'ONcells',fig13,'b','none','o')
    addLineToAxis(sparsity.controlCenterSurround(OFFcellInds),sparsity.shuffleCenterSurround(OFFcellInds),...
        'OFFcells',fig13,'r','none','o')
    downLim = min([sparsity.controlCenterSurround, sparsity.shuffleCenterSurround, 0]);
    upLim = 1.1*max([sparsity.controlCenterSurround, sparsity.shuffleCenterSurround]);
    addLineToAxis([downLim upLim],[downLim upLim],'unity',fig13,'k','--','none')
    
    
    recID = getRecordingTypeFromEpochList(currentNode.epochList);
    if (exportFigs)
        figID = ['CSNIL_controlInd_',recID];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')

        figID = ['CSNIL_controlBoth_',recID];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')

        figID = ['CSNIL_shuffInd_',recID];
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')

        figID = ['CSNIL_shuffBoth_',recID];
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_scatterFixations_',recID];
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_popDiff_',recID];
        makeAxisStruct(fig7,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_stimCtl_',recID];
        makeAxisStruct(fig8,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_stimShuff_',recID];
        makeAxisStruct(fig9,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_stimScatterCtl_',recID];
        makeAxisStruct(fig10,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_stimScatterShuff_',recID];
        makeAxisStruct(fig11,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_imageHistogram_',recID];
        makeAxisStruct(fig12,figID ,'RFSurroundFigs')
        
        figID = ['CSNIL_sparsity_',recID];
        makeAxisStruct(fig13,figID ,'RFSurroundFigs')
    end
end





%PORTED OVER THIS CODE FROM doCSLNAnalysis. Seems like it belongs here more
%natural image luminances in 2D space
% %             load('~/Documents/MATLAB/turner-package/resources/SaccadeLuminanceTrajectoryStimuli_20160919.mat')
% %             numberOfBins_em = 100^2;
% %             centerGenSignal = [];
% %             surroundGenSignal = [];
% %             allCStim = [];
% %             allSStim = [];
% %             for ss = 1:length(luminanceData)
% %                 cStim = resample(luminanceData(ss).centerTrajectory,center.sampleRate,200);
% %                 cStim = (cStim) ./ luminanceData(ss).ImageMax; %stim as presented
% % 
% %                 %convert to contrast (relative to mean) for filter convolution
% %                 imMean = (luminanceData(ss).ImageMean  ./ luminanceData(ss).ImageMax);
% %                 cStim = (cStim - imMean) / imMean;
% %                 allCStim = cat(2,allCStim,cStim);
% % 
% %                 linearPrediction = conv(cStim,center.LinearFilter);
% %                 linearPrediction = linearPrediction(1:length(cStim));
% %                 centerGenSignal = cat(2,centerGenSignal,linearPrediction);
% % 
% %                 sStim = resample(luminanceData(ss).surroundTrajectory,surround.sampleRate,200);
% % 
% %                 sStim = (sStim) ./ luminanceData(ss).ImageMax; %stim as presented
% % 
% %                 %convert to contrast (relative to mean) for filter convolution
% %                 imMean = (luminanceData(ss).ImageMean  ./ luminanceData(ss).ImageMax);
% %                 sStim = (sStim - imMean) / imMean;
% %                 allSStim = cat(2,allSStim,sStim);
% % 
% %                 linearPrediction = conv(sStim,surround.LinearFilter);
% %                 linearPrediction = linearPrediction(1:length(sStim));
% %                 surroundGenSignal = cat(2,surroundGenSignal,linearPrediction);
% %             end
% % 
% %             figure; clf; fig19=gca; %2D stimulus space with eye movements
% %             set(gcf, 'WindowStyle', 'docked')
% %             set(fig19,'XScale','linear','YScale','linear')
% %             set(0, 'DefaultAxesFontSize', 12)
% %             set(get(fig19,'XLabel'),'String','Center gen. signal')
% %             set(get(fig19,'YLabel'),'String','Surround gen. signal')
% %             set(get(fig19,'YLabel'),'String','Probability')
% %             [N,Xedges,Yedges] = histcounts2(centerGenSignal,surroundGenSignal,sqrt(numberOfBins_em),...
% %                 'Normalization','probability');
% %             Ccenters = Xedges(1:end-1) + diff(Xedges);
% %             Scenters = Yedges(1:end-1) + diff(Yedges);
% %             surf(Ccenters,Scenters,log10(N))
% % %             histogram2(centerGenSignal,surroundGenSignal,sqrt(numberOfBins_em),...
% % %                 'Normalization','probability','ShowEmptyBins','on','FaceColor','flat');
% %             colormap(hot)
% %             xlabel('Center'); ylabel('Surround'); zlabel('Probability')
% %             colorbar