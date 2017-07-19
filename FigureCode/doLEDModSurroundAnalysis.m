function doLEDModSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'figureID',[],@ischar);
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    metric = ip.Results.metric;
    figureID = ip.Results.figureID;
    
    egTraceSurroundContrast = 0.75;
    
    %Image vs Disc - e.g.
    figure; clf;
    fig2=gca;
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Image center')
    set(get(fig2,'YLabel'),'String','Disc center')
    set(gcf, 'WindowStyle', 'docked')
    
    %NLI vs surround contrast - POP
    figure; clf;
    fig3=gca;
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Surround contrast')
    set(get(fig3,'YLabel'),'String','Mean NLI')
    set(gcf, 'WindowStyle', 'docked')
    
    %NLI for no surround vs natural surround intensity - pop
    figure; clf;
    fig4=gca;
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','NLI no surround')
    set(get(fig4,'YLabel'),'String','NLI with natural surround')
    set(gcf, 'WindowStyle', 'docked')
    
    %scatter with connecting lines for +/- natural surround - pop
    figure; clf;
    fig5=gca;
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','Image center')
    set(get(fig5,'YLabel'),'String','Disc center')
    set(gcf, 'WindowStyle', 'docked')
    
% % %     %Spike diff. vs surround contrast - POP
% % %     figure; clf;
% % %     fig6=gca;
% % %     set(fig6,'XScale','linear','YScale','linear')
% % %     set(0, 'DefaultAxesFontSize', 12)
% % %     set(get(fig6,'XLabel'),'String','Surround contrast')
% % %     set(get(fig6,'YLabel'),'String','Diff (spikes)')
% % %     set(gcf, 'WindowStyle', 'docked')
    
    
    %Raster - s = 0, image
    figure; clf;
    fig7=gca;
    set(fig7,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig7,'XLabel'),'String','Time (s)')
    set(get(fig7,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    %Raster - s = 0, disc
    figure; clf;
    fig8=gca;
    set(fig8,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig8,'XLabel'),'String','Time (s)')
    set(get(fig8,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    %Raster - s = eg surround, image
    figure; clf;
    fig9=gca;
    set(fig9,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig9,'XLabel'),'String','Time (s)')
    set(get(fig9,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    %Raster - s = eg surround, disc
    figure; clf;
    fig10=gca;
    set(fig10,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig10,'XLabel'),'String','Time (s)')
    set(get(fig10,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    %Time to first spike, s = +0.9
    figure; clf;
    fig11=gca;
    set(fig11,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig11,'XLabel'),'String','Surround contrast')
    set(get(fig11,'YLabel'),'String','Time to first spike (ms)')
    set(gcf, 'WindowStyle', 'docked')
    
    %KNN discrim. pCorrect vs surround contrast
    figure; clf;
    fig12=gca;
    set(fig12,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig12,'XLabel'),'String','Surround contrast')
    set(get(fig12,'YLabel'),'String','pCorrect(Image vs. Disc)')
    set(gcf, 'WindowStyle', 'docked')
    
    %NLI vs (Ic - Is)
    figure; clf;
    fig13=gca;
    set(fig13,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig13,'XLabel'),'String','Ic - Is')
    set(get(fig13,'YLabel'),'String','NLI')
    set(gcf, 'WindowStyle', 'docked')
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(imageName)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    
    NLIvalues.paramSurrounds = [];
    RespDiff.paramSurrounds = [];
    NLIvalues.naturalSurround = [];
    resp_noSurround = [];
    resp_natSurround = [];
    timeToFirstSpike.image = [];
    timeToFirstSpike.disc = [];
    
    runningPCorrect.timingCode = [];
    runningPCorrect.countCode = [];
    
    AllIntDiff = [];
    AllNLI = [];
    AllRespDiff = [];
    
    pCorrect_natural = [];
    
    patchCt = 0;
    for pp = 1:length(populationNodes)
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        for ii = 1:cellNode.children.length %for each image
            imageNode = cellNode.children(ii);
            noPatches = imageNode.children.length;
            noSurrounds = imageNode.children(1).children.length;
            %rows: patch 1 image, patch 1 disc, patch 2 image, ...
            %cols: surrounds, [parameterized array, image patch surround]
            clear imageResponseMatrix
            imageResponseMatrix.image.mean = nan(noPatches,noSurrounds-1);
            imageResponseMatrix.image.err = nan(noPatches,noSurrounds-1);
            imageResponseMatrix.disc.mean = nan(noPatches,noSurrounds-1);
            imageResponseMatrix.disc.err = nan(noPatches,noSurrounds-1);
            pCorrect.timingCode = nan(noPatches,noSurrounds-1);
            pCorrect.countCode = nan(noPatches,noSurrounds-1);

            surroundContrastValues = nan(noPatches,noSurrounds-1);
            for ll = 1:noPatches
                patchCt = patchCt + 1;
                patchNode = imageNode.children(ll);
                parameterizedSurroundContrasts = convertJavaArrayList(...
                    patchNode.epochList.firstValue.protocolSettings('surroundContrast'));
                NLIvalues.paramSurrounds = cat(1,NLIvalues.paramSurrounds, nan(1,9));
                RespDiff.paramSurrounds = cat(1,RespDiff.paramSurrounds, nan(1,9));
                
                for ss = 1:noSurrounds
                    currentSurroundContrast = patchNode.children(ss).splitValue;

                    imageResp = getResponseAmplitudeStats(patchNode.children(ss).childBySplitValue('image').epochList,recType);
                    discResp = getResponseAmplitudeStats(patchNode.children(ss).childBySplitValue('intensity').epochList,recType);

                    newDiff = (imageResp.(metric).mean - discResp.(metric).mean);
                    newNLIvalue =  newDiff ./ ...
                           (abs(imageResp.(metric).mean) + abs(discResp.(metric).mean));
                       
                    %For spike rec: time to first spike and kNN stuff.
                    sampleRate = imageNode.epochList.firstValue.protocolSettings('sampleRate');
                    startPts = (imageNode.epochList.firstValue.protocolSettings('preTime') / 1e3) * sampleRate;
                    endPts = (imageNode.epochList.firstValue.protocolSettings('stimTime') / 1e3) * sampleRate + startPts;
                    if strcmp(recType,'extracellular')
                        imageTrace = getMeanResponseTrace(patchNode.children(ss).childBySplitValue('image').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        discTrace = getMeanResponseTrace(patchNode.children(ss).childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        tempImageTime = [];
                        for tt = 1:size(imageTrace.binary,1)
                            firstSpike = (find(imageTrace.binary(tt,startPts:endPts),1) / sampleRate) * 1e3;
                            tempImageTime = cat(2,tempImageTime,firstSpike);
                        end
                        tempDiscTime = [];
                        for tt = 1:size(discTrace.binary,1)
                            firstSpike = (find(discTrace.binary(tt,startPts:endPts),1) / sampleRate) * 1e3;
                            tempDiscTime = cat(2,tempDiscTime,firstSpike);
                        end
                        %kNN observations:
                        inputObservations = cell(1,2);
                        inputObservations{1} = imageTrace.binary(:,startPts:endPts);
                        inputObservations{2} = discTrace.binary(:,startPts:endPts);
                        NumNeighbors = 3;
                    end
            
                    if ismember(currentSurroundContrast,parameterizedSurroundContrasts) %non-natural, param surround
                       putInd = find(currentSurroundContrast==parameterizedSurroundContrasts);
                       NLIvalues.paramSurrounds(patchCt,putInd) = newNLIvalue;
                       RespDiff.paramSurrounds(patchCt,putInd) = newDiff;
                       
                       imageResponseMatrix.image.mean(ll,putInd) = imageResp.(metric).mean;
                       imageResponseMatrix.image.err(ll,putInd) = imageResp.(metric).SEM;

                       imageResponseMatrix.disc.mean(ll,putInd) = discResp.(metric).mean;
                       imageResponseMatrix.disc.err(ll,putInd) = discResp.(metric).SEM;
                       surroundContrastValues(ll,putInd) = currentSurroundContrast;
                       
                        %spike timing matters:
                        precision=40; %msec
                        precision=precision*10; %datapoints
                        q_cost=2/precision;
                        res_timing = doKNNBySpikeD(inputObservations,NumNeighbors,q_cost);
                        pCorrect.timingCode(ll,putInd) = res_timing.pCorrect;

                        %just spike count:
                        q_cost=0; %free to move spikes. i.e. spike count code
                        res_count = doKNNBySpikeD(inputObservations,NumNeighbors,q_cost);
                        pCorrect.countCode(ll,putInd) = res_count.pCorrect;

                        timeToFirstSpike.image(patchCt,putInd) = median(tempImageTime);
                        timeToFirstSpike.disc(patchCt,putInd) = median(tempDiscTime);
                    else %natural surround
                        NLIvalues.naturalSurround(patchCt) = newNLIvalue;
                        resp_natSurround(patchCt,:) = [imageResp.(metric).mean, discResp.(metric).mean];
                        
                        %just spike count:
                        q_cost=0; %free to move spikes. i.e. spike count code
                        res_count = doKNNBySpikeD(inputObservations,NumNeighbors,q_cost);
                        pCorrect_natural(patchCt,1) = res_count.pCorrect;
                       
                        %spike timing matters:
                        precision=40; %msec
                        precision=precision*10; %datapoints
                        q_cost=2/precision;
                        res_timing = doKNNBySpikeD(inputObservations,NumNeighbors,q_cost);
                        pCorrect_natural(patchCt,2) = res_timing.pCorrect;

                        timeToFirstSpike.image(patchCt,putInd) = median(tempImageTime);
                        timeToFirstSpike.disc(patchCt,putInd) = median(tempDiscTime);
                       
                    end
                    if currentSurroundContrast == 0
                       resp_noSurround(patchCt,:) = [imageResp.(metric).mean, discResp.(metric).mean];
                    end

                    if patchNode.custom.get('isExample')
                        addLineToAxis(0,0,cellInfo.cellID,fig2,'k','none','none')
                        addLineToAxis(0,0,cellInfo.cellID,fig7,'k','none','none')
                        addLineToAxis(0,0,cellInfo.cellID,fig8,'k','none','none')
                        addLineToAxis(0,0,cellInfo.cellID,fig9,'k','none','none')
                        addLineToAxis(0,0,cellInfo.cellID,fig10,'k','none','none')
                        if (currentSurroundContrast == egTraceSurroundContrast)
                            %add example traces, with target surround
                            imageTrace = getMeanResponseTrace(patchNode.children(ss).childBySplitValue('image').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                            discTrace = getMeanResponseTrace(patchNode.children(ss).childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);

                            if strcmp(recType,'extracellular') %raster plot - eg surround
                                addRastersToFigure(imageTrace.binary,fig9)
                                addRastersToFigure(discTrace.binary,fig10)
                            end

                        elseif (currentSurroundContrast == 0)
                            %add example traces, no surround
                            imageTrace = getMeanResponseTrace(patchNode.children(ss).childBySplitValue('image').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                            discTrace = getMeanResponseTrace(patchNode.children(ss).childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);

                            if strcmp(recType,'extracellular') %raster plot - no surround
                                addRastersToFigure(imageTrace.binary,fig7)
                                addRastersToFigure(discTrace.binary,fig8)
                            end %end if raster
                            
                        end %end if eg surrounds
                    end %end if example
                end %for surround


                %NLI vs (Ic - Is)
                backgroundInt = patchNode.epochList.firstValue.protocolSettings('backgroundIntensity');
                discInt = patchNode.epochList.firstValue.protocolSettings('equivalentIntensity');
                surroundInt = backgroundInt + backgroundInt * parameterizedSurroundContrasts;
                IntDiff = discInt - surroundInt;
                IntDiff = IntDiff ./ backgroundInt; %normalize by image background
                
                tempDiff = (imageResponseMatrix.image.mean(ll,:) - imageResponseMatrix.disc.mean(ll,:));
                tempNLI = tempDiff ./ ...
                    (imageResponseMatrix.image.mean(ll,:) + imageResponseMatrix.disc.mean(ll,:));
                
                
                
                AllIntDiff(patchCt,:) = IntDiff;
                AllNLI(patchCt,:) = tempNLI;
                AllRespDiff(patchCt,:) = tempDiff;
                
                
               if patchNode.custom.get('isExample') %add to example figs
                    patchLocation = str2num(patchNode.splitValue);
                    imageRes = getNaturalImagePatchFromLocation(patchLocation,imageNode.splitValue,'imageSize',[250 250]);
                    figure(30); imagesc(imageRes.images{1}); colormap(gray); axis image; axis off;
                    presentedContrasts = surroundContrastValues(ll,:);
                    colors = pmkmp(length(presentedContrasts));

                    for cc = 1:length(imageResponseMatrix.image.mean(ll,:))
                        plotMarker = 'o';
                        addLineToAxis(imageResponseMatrix.image.mean(ll,cc),...
                            imageResponseMatrix.disc.mean(ll,cc),...
                            ['mean',num2str(cc)],fig2,colors(cc,:),'none',plotMarker)
                        tempX = imageResponseMatrix.image.mean(ll,cc) + ...
                            [-imageResponseMatrix.image.err(ll,cc), imageResponseMatrix.image.err(ll,cc)];
                        tempY = [imageResponseMatrix.disc.mean(ll,cc), imageResponseMatrix.disc.mean(ll,cc)];
                        addLineToAxis(tempX, tempY,['errX',num2str(cc)],fig2,colors(cc,:),'-','none')
                        
                        tempX = [imageResponseMatrix.image.mean(ll,cc), imageResponseMatrix.image.mean(ll,cc)];
                        tempY = imageResponseMatrix.disc.mean(ll,cc) + ...
                            [-imageResponseMatrix.disc.err(ll,cc), imageResponseMatrix.disc.err(ll,cc)];
                        addLineToAxis(tempX, tempY,['errY',num2str(cc)],fig2,colors(cc,:),'-','none')
                    end
                    
                    limUp = max([imageResponseMatrix.image.mean(ll,:), imageResponseMatrix.disc.mean(ll,:)]);
                    limDown = min([imageResponseMatrix.image.mean(ll,:), imageResponseMatrix.disc.mean(ll,:), 0]);
                    addLineToAxis([limDown limUp],[limDown limUp],...
                            'unity',fig2,'k','--','none')

                end %if example plotting
            end %for patch location
        
            runningPCorrect.timingCode = cat(1,runningPCorrect.timingCode,pCorrect.timingCode);
            runningPCorrect.countCode = cat(1,runningPCorrect.countCode,pCorrect.countCode);

        end %for image

    end %for cell

    %NLI vs (Ic - Is)
    binAndPlotEquallyPopulatedBins(AllIntDiff(:),AllNLI(:), 9, fig13)
    

    
    disp('-------------------')
    disp([num2str(length(NLIvalues.naturalSurround)), ' patches'])
    disp([num2str(pp), ' cells'])
    
    %param surrounds pCorrect:
    addLineToAxis(parameterizedSurroundContrasts,mean(runningPCorrect.countCode),'count',fig12,'k','-','.')
    tempErr = std(runningPCorrect.countCode,[],1) ./ sqrt(size(runningPCorrect.countCode,1));
    addLineToAxis(parameterizedSurroundContrasts,mean(runningPCorrect.countCode) + tempErr,'counteUp',fig12,'k','--','none')
    addLineToAxis(parameterizedSurroundContrasts,mean(runningPCorrect.countCode) - tempErr,'counteDown',fig12,'k','--','none')

    addLineToAxis(parameterizedSurroundContrasts,mean(runningPCorrect.timingCode),'timing',fig12,'r','-','.')
    tempErr = std(runningPCorrect.timingCode,[],1) ./ sqrt(size(runningPCorrect.timingCode,1));
    addLineToAxis(parameterizedSurroundContrasts,mean(runningPCorrect.timingCode) + tempErr,'timingeUp',fig12,'r','--','none')
    addLineToAxis(parameterizedSurroundContrasts,mean(runningPCorrect.timingCode) - tempErr,'timingeDown',fig12,'r','--','none')


    %pCorrect nat. surround, timing vs count:
    [~, p] = ttest(pCorrect_natural(:,1),pCorrect_natural(:,2));
    disp('Nat. surround:')
    disp(mean(pCorrect_natural))
    disp(std(pCorrect_natural) ./ sqrt(patchCt))
    disp('p = ')
    disp(p);

    % pop. NLI vs. surround contrast:
    tempMean = mean(NLIvalues.paramSurrounds,1);
    tempErr = std(NLIvalues.paramSurrounds,[],1) ./ sqrt(size(NLIvalues.paramSurrounds,1));
    addLineToAxis(parameterizedSurroundContrasts,tempMean,...
        'nliMean',fig3,'k','none','o')
    addLineToAxis(parameterizedSurroundContrasts,tempMean + tempErr,...
        'nliErrUp',fig3,'k','--','none')
    addLineToAxis(parameterizedSurroundContrasts,tempMean - tempErr,...
        'nliErrDown',fig3,'k','--','none')

    % pop. response difference vs. surround contrast:
    tempMean = mean(RespDiff.paramSurrounds,1);
    tempErr = std(RespDiff.paramSurrounds,[],1) ./ sqrt(size(RespDiff.paramSurrounds,1));
    addLineToAxis(parameterizedSurroundContrasts,tempMean,...
        'diffMean',fig3,'r','none','o')
    addLineToAxis(parameterizedSurroundContrasts,tempMean + tempErr,...
        'diffErrUp',fig3,'r','--','none')
    addLineToAxis(parameterizedSurroundContrasts,tempMean - tempErr,...
        'diffErrDown',fig3,'r','--','none')

    
    
    zeroInd = parameterizedSurroundContrasts == 0;
    addLineToAxis(NLIvalues.paramSurrounds(:,zeroInd),NLIvalues.naturalSurround, 'dat',fig4,'k','none','o')
    
    tempMeanX = mean(NLIvalues.paramSurrounds(:,zeroInd));
    tempMeanY = mean(NLIvalues.naturalSurround);
    errX = std(NLIvalues.paramSurrounds(:,zeroInd)) / sqrt(length(NLIvalues.paramSurrounds(:,zeroInd)));
    errY = std(NLIvalues.naturalSurround) / sqrt(length(NLIvalues.naturalSurround));
    addLineToAxis(tempMeanX, tempMeanY, 'meanPt',fig4,'k','none','.')
    addLineToAxis(tempMeanX + [errX, -errX], [tempMeanY tempMeanY], 'errX',fig4,'k','-','none')
    addLineToAxis([tempMeanX, tempMeanX], tempMeanY + [errY, -errY], 'errY',fig4,'k','-','none')
    
    limDown = min([NLIvalues.paramSurrounds(:)', NLIvalues.naturalSurround(:)', 0]);
    limUp = 1.1*max([NLIvalues.paramSurrounds(:)', NLIvalues.naturalSurround(:)']);
    addLineToAxis([limDown limUp],[limDown limUp],'unity',fig4,'k','--','none')
    

    addLineToAxis(resp_noSurround(:,1),resp_noSurround(:,2),...
        'noSurround',fig5,colors(1,:),'none','.')
    addLineToAxis(resp_natSurround(:,1),resp_natSurround(:,2),...
        'natSurround',fig5,colors(5,:),'none','o')
    for pp = 1:patchCt
        addLineToAxis([resp_noSurround(pp,1), resp_natSurround(pp,1)],...
            [resp_noSurround(pp,2), resp_natSurround(pp,2)],...
        ['lineTo',num2str(pp)],fig5,[0.3 0.3 0.3],'-','none')
    end
    limDown = min([resp_natSurround(:)', resp_natSurround(:)', 0]);
    limUp = 1.1*max([resp_natSurround(:)', resp_natSurround(:)']);
    addLineToAxis([limUp limDown],[limUp limDown],...
        'unity',fig5,'k','--','none')
    
    %Time to first spike:
    meanIm = nanmean(timeToFirstSpike.image,1);
    errIm = nanstd(timeToFirstSpike.image,[],1) ./ sqrt(size(timeToFirstSpike.image,1));
    addLineToAxis(parameterizedSurroundContrasts,meanIm,'imageMean',fig11,'k','none','o')
    addLineToAxis(parameterizedSurroundContrasts,meanIm+errIm,'imageErrUp',fig11,'k','--','none')
    addLineToAxis(parameterizedSurroundContrasts,meanIm-errIm,'imageErrDown',fig11,'k','--','none')
    
    meanDisc = nanmean(timeToFirstSpike.disc,1);
    errDisc = nanstd(timeToFirstSpike.disc,[],1) ./ sqrt(size(timeToFirstSpike.disc,1));
    addLineToAxis(parameterizedSurroundContrasts,meanDisc,'discMean',fig11,'g','none','o')
    addLineToAxis(parameterizedSurroundContrasts,meanDisc+errDisc,'discErrUp',fig11,'g','--','none')
    addLineToAxis(parameterizedSurroundContrasts,meanDisc-errDisc,'discErrDown',fig11,'g','--','none')
    
    if ~isempty(figureID)
        makeAxisStruct(fig2,['LEDmodSeg_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig3,['LEDmodS_NLIvsS_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig4,['LEDmodS_NLInat_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig5,['LEDmodS_natScatter_',figureID] ,'RFSurroundFigs')
%         makeAxisStruct(fig6,['LEDmodS_DvsS_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig13,['LEDmodS_IntDiff_',figureID] ,'RFSurroundFigs')
        if strcmp(recType,'extracellular') %raster plots
            figID = 'LEDmodS_rast_imAlone';
            makeAxisStruct(fig7,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_rast_discAlone';
            makeAxisStruct(fig8,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_rast_imSurr';
            makeAxisStruct(fig9,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_rast_discSurr';
            makeAxisStruct(fig10,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_firstSpike';
            makeAxisStruct(fig11,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_KNN_pCorr';
            makeAxisStruct(fig12,figID ,'RFSurroundFigs')
        end
    end
end