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
    
    %NLI vs surround contrast - POP
    figure; clf;
    fig3=gca;
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Surround contrast')
    set(get(fig3,'YLabel'),'String','Mean NLI')
    
    %NLI for no surround vs natural surround intensity - pop
    figure; clf;
    fig4=gca;
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','NLI no surround')
    set(get(fig4,'YLabel'),'String','NLI with natural surround')
    
    %scatter with connecting lines for +/- natural surround - pop
    figure; clf;
    fig5=gca;
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','Image center')
    set(get(fig5,'YLabel'),'String','Disc center')
    
    
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
    NLIvalues.naturalSurround = [];
    resp_noSurround = [];
    resp_natSurround = [];
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
            surroundContrastValues = nan(noPatches,noSurrounds-1);
            for ll = 1:noPatches
               patchCt = patchCt + 1;
               patchNode = imageNode.children(ll);
               parameterizedSurroundContrasts = convertJavaArrayList(...
                        patchNode.epochList.firstValue.protocolSettings('surroundContrast'));
               NLIvalues.paramSurrounds = cat(1,NLIvalues.paramSurrounds, nan(1,9));  
                for ss = 1:noSurrounds
                    currentSurroundContrast = patchNode.children(ss).splitValue;

                    imageResp = getResponseAmplitudeStats(patchNode.children(ss).childBySplitValue('image').epochList,recType);
                    discResp = getResponseAmplitudeStats(patchNode.children(ss).childBySplitValue('intensity').epochList,recType);

                    newNLIvalue = (imageResp.(metric).mean - discResp.(metric).mean) ./ ...
                           (abs(imageResp.(metric).mean) + abs(discResp.(metric).mean));
                                          
                    if ismember(currentSurroundContrast,parameterizedSurroundContrasts)
                       putInd = find(currentSurroundContrast==parameterizedSurroundContrasts);
                       NLIvalues.paramSurrounds(patchCt,putInd) = newNLIvalue;

                       imageResponseMatrix.image.mean(ll,putInd) = imageResp.(metric).mean;
                       imageResponseMatrix.image.err(ll,putInd) = imageResp.(metric).SEM;

                       imageResponseMatrix.disc.mean(ll,putInd) = discResp.(metric).mean;
                       imageResponseMatrix.disc.err(ll,putInd) = discResp.(metric).SEM;
                       surroundContrastValues(ll,putInd) = currentSurroundContrast;
                    else
                       NLIvalues.naturalSurround(patchCt) = newNLIvalue;
                       resp_natSurround(patchCt,:) = [imageResp.(metric).mean, discResp.(metric).mean];
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
        
        end %for image

    end %for cell
    disp('-------------------')
    disp([num2str(length(NLIvalues.naturalSurround)), ' patches'])
    disp([num2str(pp), ' cells'])
    for ss = 1:length(parameterizedSurroundContrasts)
        currentContrast = parameterizedSurroundContrasts(ss);
        meanPts = mean(NLIvalues.paramSurrounds(:,ss));
        errPts = std(NLIvalues.paramSurrounds(:,ss)) ./ sqrt(length(NLIvalues.paramSurrounds(:,ss)));
        addLineToAxis(currentContrast,meanPts,...
            ['mean',num2str(ss)],fig3,'k','none','o')
        addLineToAxis([currentContrast, currentContrast],[meanPts - errPts, meanPts + errPts],...
            ['err',num2str(ss)],fig3,'k','-','none')
    end
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
    
    if ~isempty(figureID)
        makeAxisStruct(fig2,['LEDmodSeg_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig3,['LEDmodS_NLIvsS_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig4,['LEDmodS_NLInat_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig5,['LEDmodS_natScatter_',figureID] ,'RFSurroundFigs')
        
        if strcmp(recType,'extracellular') %raster plots
            figID = 'LEDmodS_rast_imAlone';
            makeAxisStruct(fig7,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_rast_discAlone';
            makeAxisStruct(fig8,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_rast_imSurr';
            makeAxisStruct(fig9,figID ,'RFSurroundFigs')
            
            figID = 'LEDmodS_rast_discSurr';
            makeAxisStruct(fig10,figID ,'RFSurroundFigs')
        end
    end
end