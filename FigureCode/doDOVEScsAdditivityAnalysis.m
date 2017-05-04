function doDOVEScsAdditivityAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    
    figColors = pmkmp(8);
    figDir = '~/Documents/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs
    
    
    figure; clf; fig1=gca; initFig(fig1,'Time (s)','Response') % Response traces, c and s, spikes
    figure; clf; fig2=gca; initFig(fig2,'Time (s)','Response') % Response traces, linear sum and cs, spikes
    
    figure; clf; fig3=gca; initFig(fig3,'Time (s)','Response') % Response traces, c and s, exc
    figure; clf; fig4=gca; initFig(fig4,'Time (s)','Response') % Response traces, linear sum and cs, exc
    
    figure; clf; fig5=gca; initFig(fig5,'Time (s)','Location (arcmin)') % Eye movement positions over time
    
    figure; clf; fig6=gca; initFig(fig6,'Diff. Exc. (pA)','Diff. Spikes') % Spikes vs. Exc difference
    
    figure; clf; fig7=gca; initFig(fig7,'Measured (spikes)','Lin. sum (spikes)') % Pop: spikes lin sum vs measured
    figure; clf; fig8=gca; initFig(fig8,'Measured (pC)','Lin. sum (pC)') % Pop: exc lin sum vs measured
    figure; clf; fig9=gca; initFig(fig9,'Cell type','Correlation coefficient') % Pop: diff correlations
    
    figure; clf; fig11=gca; initFig(fig11,'Time (s)','Trial') % RASTER - center alone (spike rec only)
    figure; clf; fig12=gca; initFig(fig12,'Time (s)','Trial') % RASTER - surround alone (spike rec only)
    figure; clf; fig13=gca; initFig(fig13,'Time (s)','Trial') % RASTER - center + surround (spike rec only)
    
    populationNodes = {}; %pop. node is a stimulus
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                '@(list)splitOnRecKeyword(list)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    diffCorrValues = [];
    meanVals.spikes.linsum = [];
    meanVals.spikes.measured = [];

    meanVals.exc.linsum = [];
    meanVals.exc.measured = [];
    ONcellInds = [];
    OFFcellInds = [];
    for pp = 1:length(populationNodes) %for cell in population
        currentNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(currentNode.epochList);
        if strcmp(cellInfo.cellType,'ONparasol')
            ONcellInds = cat(2,ONcellInds,pp);
        elseif strcmp(cellInfo.cellType,'OFFparasol')
            OFFcellInds = cat(2,OFFcellInds,pp);
        end
        %exc:
        excNode = currentNode.childBySplitValue('exc');
        recType = getRecordingTypeFromEpochList(excNode.epochList);
        exc_centerRes = getMeanResponseTrace(excNode.childBySplitValue('Center').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        exc_surroundRes = getMeanResponseTrace(excNode.childBySplitValue('Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        exc_centerSurroundRes = getMeanResponseTrace(excNode.childBySplitValue('Center-Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        %spikes:
        spikesNode = currentNode.childBySplitValue('extracellular');
        recType = getRecordingTypeFromEpochList(spikesNode.epochList);
        spikes_centerRes = getMeanResponseTrace(spikesNode.childBySplitValue('Center').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        spikes_surroundRes = getMeanResponseTrace(spikesNode.childBySplitValue('Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        spikes_centerSurroundRes = getMeanResponseTrace(spikesNode.childBySplitValue('Center-Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        
        %timing stuff:
        timingEpoch = currentNode.epochList.firstValue;
        frameRate = timingEpoch.protocolSettings('background:Microdisplay Stage@localhost:monitorRefreshRate');
        if isempty(frameRate)
            frameRate = timingEpoch.protocolSettings('background:Microdisplay_Stage@localhost:monitorRefreshRate');
        end
        FMdata = (riekesuite.getResponseVector(timingEpoch,'Frame Monitor'))';
        [measuredFrameTimes, ~] = getFrameTiming(FMdata,0);

        sampleRate = currentNode.epochList.firstValue.protocolSettings('sampleRate');
        preTime = currentNode.epochList.firstValue.protocolSettings('preTime') / 1e3; %sec
        stimTime = currentNode.epochList.firstValue.protocolSettings('stimTime') / 1e3; %sec
        preFrames = preTime .* frameRate;
        firstStimFrameFlip = measuredFrameTimes(preFrames + 1);
        
        %Stimulus stuff. To get fixation periods...
        stimSet = currentNode.epochList.firstValue.protocolSettings('currentStimSet');
        imageSet = 'NaturalImages';
        resourcesDir = '~/Documents/MATLAB/turner-package/resources/';
        load([resourcesDir, stimSet]);
        eyeX = FEMdata(currentNode.splitValue).eyeX;
        eyeY = FEMdata(currentNode.splitValue).eyeY;
        dx = abs(diff(eyeX));
        dy = abs(diff(eyeY));
        instVel = sqrt(dx.^2 + dy.^2); %instantaneous velocity
        
        [peaks,peakInd] = getPeaks(instVel,1);
        threshInd = peaks > 10; %only changes greater than 10 amin / unit time
        peakInd = peakInd(threshInd);
        
        %10 sample points at 200 Hz -> 50 msec apart. Filters out quick
        %little "compound" saccadic movements
        refractoryPeriod = 10;
        filtInds = peakInd(1);
        ct = 1;
        for ii = 2:length(peakInd)
            if peakInd(ii) - filtInds(ct) < refractoryPeriod
            else
                ct = ct + 1;
                filtInds = cat(1,filtInds,peakInd(ii));
            end
        end
        filtInds = [0 filtInds']; %start at 1 and shift back due to diff above
        filtInds(end + 1) = length(eyeX); %end of last fixation
        
        saccadeTimes = round((filtInds ./ 200) .* sampleRate + firstStimFrameFlip); %datapoints
        fixationResponses = nan(6,length(saccadeTimes)-1); %rows go c,s,cs (spikes) then c,s,cs (exc)
        chargeMult = -1;
        for ss = 1:(length(saccadeTimes)-1)
            tempStart = saccadeTimes(ss);
            tempEnd = saccadeTimes(ss+1);
            
            fixationResponses(1,ss) = mean(sum(spikes_centerRes.binary(:,tempStart:tempEnd),2),1); %spike count
            fixationResponses(2,ss) = mean(sum(spikes_surroundRes.binary(:,tempStart:tempEnd),2),1); %spike count
            fixationResponses(3,ss) = mean(sum(spikes_centerSurroundRes.binary(:,tempStart:tempEnd),2),1); %spike count
            
            fixationResponses(4,ss) = chargeMult * trapz(exc_centerRes.mean(tempStart:tempEnd)) / sampleRate; %pC
            fixationResponses(5,ss) = chargeMult * trapz(exc_surroundRes.mean(tempStart:tempEnd)) / sampleRate; %pC
            fixationResponses(6,ss) = chargeMult * trapz(exc_centerSurroundRes.mean(tempStart:tempEnd)) / sampleRate; %pC
        end
        LinSum_spikes = fixationResponses(1,:) + fixationResponses(2,:);
        Meas_spikes = fixationResponses(3,:);
        diff_spikes = LinSum_spikes - Meas_spikes;
        meanVals.spikes.linsum(pp) = mean(LinSum_spikes);
        meanVals.spikes.measured(pp) = mean(Meas_spikes);

        LinSum_exc = fixationResponses(4,:) + fixationResponses(5,:);
        Meas_exc = fixationResponses(6,:);
        diff_exc = LinSum_exc - Meas_exc;
        meanVals.exc.linsum(pp) = mean(LinSum_exc);
        meanVals.exc.measured(pp) = mean(Meas_exc);

        newCorr = corr(diff_exc',diff_spikes');
        diffCorrValues(pp) = newCorr;

        if currentNode.custom.get('isExample')
            %difference correlation
            addLineToAxis(diff_exc,diff_spikes,...
            ['d',num2str(pp),'_', num2str(currentNode.splitValue)],fig6,'k','none','o')
            addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(currentNode.splitValue)],fig6,'k','none','none')
            p = polyfit(diff_exc,diff_spikes,1);
            xFit = [min(diff_exc), max(diff_exc)];
            yFit = polyval(p,xFit);
            addLineToAxis(xFit,yFit,['fit',num2str(pp),'_', num2str(currentNode.splitValue)],fig6,'k','-','none')
            %corr stat for example cell:
            [rho, pval] = corr(diff_exc',diff_spikes');
            disp('Diff corr for e.g. cell')
            disp([rho, pval])
            
            
            %spikes traces:
            addLineToAxis(spikes_centerRes.timeVector,spikes_centerRes.mean,...
                'center',fig1,figColors(1,:),'-','none')
            addLineToAxis(spikes_surroundRes.timeVector,spikes_surroundRes.mean,...
                'surround',fig1,figColors(4,:),'-','none')
            addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(currentNode.splitValue)],fig1,'k','none','none')

            addLineToAxis(spikes_centerRes.timeVector,spikes_centerRes.mean + spikes_surroundRes.mean,...
                'linSum',fig2,[0.7 0.7 0.7],'-','none')
            addLineToAxis(spikes_centerSurroundRes.timeVector,spikes_centerSurroundRes.mean,...
                'measuredCS',fig2,'k','-','none')
            addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(currentNode.splitValue)],fig2,'k','none','none')
            if strcmp(recType,'extracellular') %do raster figs
                addRastersToFigure(spikes_centerRes.binary,fig11)                    
                addRastersToFigure(spikes_surroundRes.binary,fig12)                    
                addRastersToFigure(spikes_centerSurroundRes.binary,fig13)
            end
            
            %exc traces:
            addLineToAxis(exc_centerRes.timeVector,exc_centerRes.mean,...
                'center',fig3,figColors(1,:),'-','none')
            addLineToAxis(exc_surroundRes.timeVector,exc_surroundRes.mean,...
                'surround',fig3,figColors(4,:),'-','none')
            addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(currentNode.splitValue)],fig3,'k','none','none')

            addLineToAxis(exc_centerRes.timeVector,exc_centerRes.mean + exc_surroundRes.mean,...
                'linSum',fig4,[0.7 0.7 0.7],'-','none')
            addLineToAxis(exc_centerSurroundRes.timeVector,exc_centerSurroundRes.mean,...
                'measuredCS',fig4,'k','-','none')
            addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(currentNode.splitValue)],fig4,'k','none','none')

            %natural image with overlayed eye movements
            imageName = FEMdata(currentNode.splitValue).ImageName;
            fileId=fopen([resourcesDir, imageSet, '/', imageName],'rb','ieee-be');
            img = fread(fileId, [1536,1024], 'uint16');
            img = double(img);
            img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
            fh = figure(21); clf;
            imagesc(img'); colormap(gray); axis image; axis off; hold on;
            plot(FEMdata(currentNode.splitValue).eyeX,FEMdata(currentNode.splitValue).eyeY,'r')
            brighten(0.6) %brighten for display purposes
            drawnow;
            figID = 'egNatImage_DOVES';
            print(fh,[figDir,figID],'-depsc')

            %image patch
            patchSize = 240; %pixels
            centerX = FEMdata(currentNode.splitValue).eyeX(end);
            centerY = FEMdata(currentNode.splitValue).eyeY(end);
            imagePatch = img(round(centerX - patchSize/2 + 1) : round(centerX + patchSize/2),...
                round(centerY - patchSize/2 + 1) : round(centerY + patchSize/2));
            fh = figure(22); clf;
            imagesc(imagePatch'); colormap(gray); axis image; axis off;
            brighten(0.6) %brighten for display purposes
            drawnow;
            figID = 'egNatImagePatch_DOVES';
            print(fh,[figDir,figID],'-depsc')

            %eye position over time
            totalOnPoints = stimTime * 200;
            movementPoints = length(FEMdata(currentNode.splitValue).eyeX); %at 200 Hz
            xTrace = FEMdata(currentNode.splitValue).eyeX;
            yTrace = FEMdata(currentNode.splitValue).eyeY;
            if movementPoints < totalOnPoints
               xTrace((movementPoints + 1):totalOnPoints) = FEMdata(currentNode.splitValue).eyeX(end);
               yTrace((movementPoints + 1):totalOnPoints) = FEMdata(currentNode.splitValue).eyeY(end);
            end

            onTimeVec = preTime + (1:totalOnPoints) ./ 200;
            addLineToAxis(onTimeVec,xTrace,...
                'xTrace',fig5,figColors(2,:),'-','none')
            addLineToAxis(onTimeVec,yTrace,...
                'yTrace',fig5,figColors(6,:),'-','none')
            addLineToAxis([onTimeVec(1) onTimeVec(1)],[0 1200],...
                'imageON',fig5,[0.7 0.7 0.7],'--','none')
            addLineToAxis([onTimeVec(end) onTimeVec(end)],[0 1200],...
                'imageOFF',fig5,[0.7 0.7 0.7],'--','none')
        end %end if example
                    
    end %end for population
    
    %spikes pop:
    %   Off...
    addLineToAxis(meanVals.spikes.measured(OFFcellInds),meanVals.spikes.linsum(OFFcellInds),'OffData',fig7,'r','none','o')
    meanOFF_x = mean(meanVals.spikes.measured(OFFcellInds));
    meanOFF_y = mean(meanVals.spikes.linsum(OFFcellInds));
    errOFF_x = std(meanVals.spikes.measured(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    errOFF_y = std(meanVals.spikes.linsum(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    addLineToAxis(meanOFF_x,meanOFF_y,'OffMean',fig7,'r','none','.')
    addLineToAxis(meanOFF_x + [errOFF_x, -errOFF_x],[meanOFF_y meanOFF_y],'OffErrX',fig7,'r','-','none')
    addLineToAxis([meanOFF_x, meanOFF_x],meanOFF_y + [errOFF_y, -errOFF_y],'OffErrY',fig7,'r','-','none')
    
    %   On...
    addLineToAxis(meanVals.spikes.measured(ONcellInds),meanVals.spikes.linsum(ONcellInds),'OnData',fig7,'b','none','o')
    meanON_x = mean(meanVals.spikes.measured(ONcellInds));
    meanON_y = mean(meanVals.spikes.linsum(ONcellInds));
    errON_x = std(meanVals.spikes.measured(ONcellInds)) ./ sqrt(length(ONcellInds));
    errON_y = std(meanVals.spikes.linsum(ONcellInds)) ./ sqrt(length(ONcellInds));
    addLineToAxis(meanON_x,meanON_y,'OnMean',fig7,'b','none','.')
    addLineToAxis(meanON_x + [errON_x, -errON_x],[meanON_y meanON_y],'OnErrX',fig7,'b','-','none')
    addLineToAxis([meanON_x, meanON_x],meanON_y + [errON_y, -errON_y],'OnErrY',fig7,'b','-','none')
    
    %stat test. ON and OFF combined (probaby split later)
    [~,p] = ttest(meanVals.spikes.measured, meanVals.spikes.linsum);
    disp('Meas. vs LinSum, spikes, ON AND OFF:')
    disp(['p = ',num2str(p)])

    upLim = max([meanVals.spikes.measured,meanVals.spikes.linsum]);
    addLineToAxis([0 upLim],[0 upLim],'unity',fig7,'k','--','none')
    
    %exc pop:
    %   Off...
    addLineToAxis(meanVals.exc.measured(OFFcellInds),meanVals.exc.linsum(OFFcellInds),'OffData',fig8,'r','none','o')
    meanOFF_x = mean(meanVals.exc.measured(OFFcellInds));
    meanOFF_y = mean(meanVals.exc.linsum(OFFcellInds));
    errOFF_x = std(meanVals.exc.measured(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    errOFF_y = std(meanVals.exc.linsum(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    addLineToAxis(meanOFF_x,meanOFF_y,'OffMean',fig8,'r','none','.')
    addLineToAxis(meanOFF_x + [errOFF_x, -errOFF_x],[meanOFF_y meanOFF_y],'OffErrX',fig8,'r','-','none')
    addLineToAxis([meanOFF_x, meanOFF_x],meanOFF_y + [errOFF_y, -errOFF_y],'OffErrY',fig8,'r','-','none')
    
    %   On...
    addLineToAxis(meanVals.exc.measured(ONcellInds),meanVals.exc.linsum(ONcellInds),'OnData',fig8,'b','none','o')
    meanON_x = mean(meanVals.exc.measured(ONcellInds));
    meanON_y = mean(meanVals.exc.linsum(ONcellInds));
    errON_x = std(meanVals.exc.measured(ONcellInds)) ./ sqrt(length(ONcellInds));
    errON_y = std(meanVals.exc.linsum(ONcellInds)) ./ sqrt(length(ONcellInds));
    addLineToAxis(meanON_x,meanON_y,'OnMean',fig8,'b','none','.')
    addLineToAxis(meanON_x + [errON_x, -errON_x],[meanON_y meanON_y],'OnErrX',fig8,'b','-','none')
    addLineToAxis([meanON_x, meanON_x],meanON_y + [errON_y, -errON_y],'OnErrY',fig8,'b','-','none')
    
    %stats:
    [~,p] = ttest(meanVals.exc.measured, meanVals.exc.linsum);
    disp('Meas. vs LinSum, exc, ON AND OFF:')
    disp(['p = ',num2str(p)])
    
    upLim = max([meanVals.exc.measured,meanVals.exc.linsum]);
    addLineToAxis([0 upLim],[0 upLim],'unity',fig8,'k','--','none')
    
    %corr pop:
    addLineToAxis(ones(1,length(OFFcellInds)),diffCorrValues(OFFcellInds),'OffData',fig9,'r','none','o')
    meanOFF = mean(diffCorrValues(OFFcellInds));
    errOFF = std(diffCorrValues(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    addLineToAxis(1.2,meanOFF,'meanOff',fig9,'r','none','.')
    addLineToAxis([1.2, 1.2],meanOFF + [errOFF -errOFF],'errOff',fig9,'r','-','none')
        
    addLineToAxis(2.*ones(1,length(ONcellInds)),diffCorrValues(ONcellInds),'OnData',fig9,'b','none','o')
    meanON = mean(diffCorrValues(ONcellInds));
    errON = std(diffCorrValues(ONcellInds)) ./ sqrt(length(ONcellInds));
    addLineToAxis(1.8,meanON,'meanOn',fig9,'b','none','.')
    addLineToAxis([1.8, 1.8],meanON + [errON -errON],'errOn',fig9,'b','-','none')
    
    if (exportFigs)
        figID = ['DOVEScs_ind_extracellular'];
        makeAxisStruct(fig1,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_CS_extracellular'];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_ind_exc'];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_CS_exc'];
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_em'];
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_egCorr'];
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_spDiff'];
        makeAxisStruct(fig7,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_excDiff'];
        makeAxisStruct(fig8,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_popCorr'];
        makeAxisStruct(fig9,figID ,'RFSurroundFigs')
        
        figID = 'DOVEScs_rasterC';
        makeAxisStruct(fig11,figID ,'RFSurroundFigs')

        figID = 'DOVEScs_rasterS';
        makeAxisStruct(fig12,figID ,'RFSurroundFigs')

        figID = 'DOVEScs_rasterCS';
        makeAxisStruct(fig13,figID ,'RFSurroundFigs')

    end
    
end