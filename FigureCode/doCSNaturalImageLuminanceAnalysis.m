function doCSNaturalImageLuminanceAnalysis(node,varargin)
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
    
    figure; clf; fig2=gca; %eg Mean trace: C, S, indep, no shuffle
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Time (s)')
    set(get(fig2,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig3=gca; %eg Mean trace: CS together and linear sum, no shuffle
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Time (s)')
    set(get(fig3,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig4=gca; %eg Mean trace: C, S, indep, with shuffle
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','Time (s)')
    set(get(fig4,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig5=gca; %eg Mean trace: CS together and linear sum, with shuffle
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','Time (s)')
    set(get(fig5,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig6=gca; %eg scatter plot of R(C) + R(S) vs R(C+S)
    set(fig6,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig6,'XLabel'),'String','R(C+S)')
    set(get(fig6,'YLabel'),'String','R(C) + R(S)')
    set(gcf, 'WindowStyle', 'docked')

    figure; clf; fig7=gca; %population. Difference (linear sum - measured) for control vs shuffled. ON and OFF cells
    set(fig7,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig7,'XLabel'),'String','Mean difference (control c,s)')
    set(get(fig7,'YLabel'),'String','Mean difference (shuffled c,s)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig8=gca; %eg. C and S stims for control
    set(fig8,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig8,'XLabel'),'String','Time(s)')
    set(get(fig8,'YLabel'),'String','Intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig9=gca; %eg. C and S stims for shuffled
    set(fig9,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig9,'XLabel'),'String','Time(s)')
    set(get(fig9,'YLabel'),'String','Intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig10=gca; %eg. C and S stim scatter, control
    set(fig10,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig10,'XLabel'),'String','Center intensity')
    set(get(fig10,'YLabel'),'String','Surround intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig11=gca; %eg. C and S stim scatter, shuffled
    set(fig11,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig11,'XLabel'),'String','Center intensity')
    set(get(fig11,'YLabel'),'String','Surround intensity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig12=gca; %eg. Nat image intensity histogram
    set(fig12,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig12,'XLabel'),'String','Intensity')
    set(get(fig12,'YLabel'),'String','Probability')
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
        elseif strcmp(cellInfo.cellType,'OFFparasol')
            OFFcellInds = cat(2,OFFcellInds,pp);
        end
        
        ImageLuminanceNode = populationNodes{pp}.childBySplitValue('CSNaturalImageLuminance');
        filterNoiseNode = populationNodes{pp}.childBySplitValue('CenterSurroundNoise').children(1).children(1);
        
% % % % % % % % GET LINEAR FILTERS % % % % % % % % % % % % % % % % 
        % center:
        center = getLinearFilterAndPrediction(filterNoiseNode.childBySplitValue('Center').epochList,recType,...
            'seedName','centerNoiseSeed');
        % surround:
        surround = getLinearFilterAndPrediction(filterNoiseNode.childBySplitValue('Surround').epochList,recType,...
            'seedName','surroundNoiseSeed');
        
% % % % % % % % DO ADDITIVITY ANALYSIS % % % % % % % % % % % % % % % %
        fixationResponses = nan(6,1); %rows are controlCenter, controlSurround,...
        fixationCount = 0;
        for imageIndex = 1:ImageLuminanceNode.children.length
            currentNode = ImageLuminanceNode.children(imageIndex);
            for ss = 1:2
                if currentNode.children(ss).splitValue == 0 %no shuffle CS
                    controlCSNode = currentNode.children(ss);
                elseif currentNode.children(ss).splitValue == 1 %shuffle CS
                    shuffleCSNode = currentNode.children(ss);
                end
            end
            controlCenter = getMeanResponseTrace(controlCSNode.childBySplitValue('Center').epochList,recType);
            controlSurround = getMeanResponseTrace(controlCSNode.childBySplitValue('Surround').epochList,recType);
            controlCenterSurround = getMeanResponseTrace(controlCSNode.childBySplitValue('Center-Surround').epochList,recType);
            
            shuffleCenter = getMeanResponseTrace(shuffleCSNode.childBySplitValue('Center').epochList,recType);
            shuffleSurround = getMeanResponseTrace(shuffleCSNode.childBySplitValue('Surround').epochList,recType);
            shuffleCenterSurround = getMeanResponseTrace(shuffleCSNode.childBySplitValue('Center-Surround').epochList,recType);
            
            % get fixation responses for this image
            sampleRate = currentNode.epochList.firstValue.protocolSettings('sampleRate');
            fixationDuration = (currentNode.epochList.firstValue.protocolSettings('fixationDuration') / 1e3) *sampleRate; %data points
            startPoint = (currentNode.epochList.firstValue.protocolSettings('preTime') / 1e3) * sampleRate; %data points
            noFixations = currentNode.epochList.firstValue.protocolSettings('stimTime') / ...
                currentNode.epochList.firstValue.protocolSettings('fixationDuration');
            for ff = 1:noFixations
                fixationCount = fixationCount + 1;
                tempStart = (ff-1)*fixationDuration + startPoint + 1;
                tempEnd = ff*fixationDuration + startPoint;
                fixationResponses(1,fixationCount) = max(controlCenter.mean(tempStart:tempEnd)); %nS
                fixationResponses(2,fixationCount) = max(controlSurround.mean(tempStart:tempEnd));
                fixationResponses(3,fixationCount) = max(controlCenterSurround.mean(tempStart:tempEnd));
                
                fixationResponses(4,fixationCount) = max(shuffleCenter.mean(tempStart:tempEnd));
                fixationResponses(5,fixationCount) = max(shuffleSurround.mean(tempStart:tempEnd));
                fixationResponses(6,fixationCount) = max(shuffleCenterSurround.mean(tempStart:tempEnd));
                
                
% %                 fixationResponses(1,fixationCount) = trapz(controlCenter.mean(tempStart:tempEnd)) / sampleRate; %nS*s
% %                 fixationResponses(2,fixationCount) = trapz(controlSurround.mean(tempStart:tempEnd)) / sampleRate;
% %                 fixationResponses(3,fixationCount) = trapz(controlCenterSurround.mean(tempStart:tempEnd)) / sampleRate;
% %                 
% %                 fixationResponses(4,fixationCount) = trapz(shuffleCenter.mean(tempStart:tempEnd)) / sampleRate;
% %                 fixationResponses(5,fixationCount) = trapz(shuffleSurround.mean(tempStart:tempEnd)) / sampleRate;
% %                 fixationResponses(6,fixationCount) = trapz(shuffleCenterSurround.mean(tempStart:tempEnd)) / sampleRate;
            end

            if currentNode.custom.get('isExample')
                %filters:
                addLineToAxis([center.filterTimeVector],[center.LinearFilter],...
                    'center',fig1,figColors(1,:),'-','none')
                addLineToAxis([surround.filterTimeVector],[surround.LinearFilter],...
                    'surround',fig1,figColors(4,:),'-','none')
                addLineToAxis(0,0,cellInfo.cellID,fig1,'k','none','none')
                
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
                
                % scatter plot with peak conductance per fixation
                % just for this example image
                egIndsToPull = (size(fixationResponses,2) - noFixations + 1):size(fixationResponses,2);
                addLineToAxis(fixationResponses(3,egIndsToPull),...
                    fixationResponses(1,egIndsToPull)+fixationResponses(2,egIndsToPull),...
                    'control',fig6,'k','none','o')
                addLineToAxis(fixationResponses(6,egIndsToPull),...
                    fixationResponses(4,egIndsToPull)+fixationResponses(5,egIndsToPull),...
                    'shuff',fig6,'r','none','o')
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
                pad = ones(1,startPoint) .* controlCSNode.epochList.firstValue.protocolSettings('backgroundIntensity');
                controlCStim_fix = convertJavaArrayList(controlCSNode.epochList.firstValue.protocolSettings('CenterIntensity'));
                controlSStim_fix = convertJavaArrayList(controlCSNode.epochList.firstValue.protocolSettings('SurroundIntensity'));
                controlCStim = [pad, kron(controlCStim_fix,ones(1,fixationDuration)), pad];
                controlSStim = [pad, kron(controlSStim_fix,ones(1,fixationDuration)), pad];

                shuffleCStim_fix = convertJavaArrayList(shuffleCSNode.epochList.firstValue.protocolSettings('CenterIntensity'));
                shuffleSStim_fix = convertJavaArrayList(shuffleCSNode.epochList.firstValue.protocolSettings('SurroundIntensity'));
                shuffleCStim = [pad, kron(shuffleCStim_fix,ones(1,fixationDuration)), pad];
                shuffleSStim = [pad, kron(shuffleSStim_fix,ones(1,fixationDuration)), pad];
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
                resourcesDir = '~/Documents/MATLAB/turner-package/resources/';
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
        meanDiff.control(pp) = mean((fixationResponses(1,:) + fixationResponses(2,:)) - fixationResponses(3,:));
        meanDiff.shuffle(pp) = mean((fixationResponses(4,:) + fixationResponses(5,:)) - fixationResponses(6,:));
    end
    
    addLineToAxis(meanDiff.control(ONcellInds),meanDiff.shuffle(ONcellInds),...
        'ONcells',fig7,'b','none','o')
    addLineToAxis(meanDiff.control(OFFcellInds),meanDiff.shuffle(OFFcellInds),...
        'OFFcells',fig7,'r','none','o')
    downLim = min([meanDiff.control, meanDiff.shuffle, 0]);
    upLim = 1.1*max([meanDiff.control, meanDiff.shuffle]);
    addLineToAxis([downLim upLim],[downLim upLim],'unity',fig7,'k','--','none')
    
    
    recID = getRecordingTypeFromEpochList(currentNode.epochList);
    if (exportFigs)
        figID = ['CSNIL_filters_',recID];
        makeAxisStruct(fig1,figID ,'RFSurroundFigs')

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