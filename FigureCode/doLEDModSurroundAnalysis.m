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
    
    figure; clf; fig2=gca; initFig(fig2,'Image center','Disc center') %Image vs Disc - e.g.
    
    figure; clf; fig3=gca; initFig(fig3,'Surround contrast','Mean NLI') %NLI vs surround contrast - POP
    figure; clf; fig4=gca; initFig(fig4,'NLI no surround','NLI with natural surround') %NLI for no surround vs natural surround intensity - pop
    figure; clf; fig6=gca; initFig(fig6,'Surround contrast','Diff (spikes)') %Resp diff vs. surround contrast - POP
    figure; clf; fig13=gca; initFig(fig13,'Ic - Is','NLI') %NLI vs (Ic - Is)
    
    figure; clf; fig7=gca; initFig(fig7,'Time (s)','Trial') %Raster - s = 0, image
    figure; clf; fig8=gca; initFig(fig8,'Time (s)','Trial') %Raster - s = 0, disc
    figure; clf; fig9=gca; initFig(fig9,'Time (s)','Trial') %Raster - s = eg surround, image
    figure; clf; fig10=gca; initFig(fig10,'Time (s)','Trial') %Raster - s = eg surround, disc

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

    AllIntDiff = [];
    AllNLI = [];
    AllRespDiff = [];
    patchMean = [];

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

                    if ismember(currentSurroundContrast,parameterizedSurroundContrasts) %non-natural, param surround
                       putInd = find(currentSurroundContrast==parameterizedSurroundContrasts);
                       NLIvalues.paramSurrounds(patchCt,putInd) = newNLIvalue;
                       RespDiff.paramSurrounds(patchCt,putInd) = newDiff;
                       
                       imageResponseMatrix.image.mean(ll,putInd) = imageResp.(metric).mean;
                       imageResponseMatrix.image.err(ll,putInd) = imageResp.(metric).SEM;

                       imageResponseMatrix.disc.mean(ll,putInd) = discResp.(metric).mean;
                       imageResponseMatrix.disc.err(ll,putInd) = discResp.(metric).SEM;
                       surroundContrastValues(ll,putInd) = currentSurroundContrast;
                    else %natural surround
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

                %NLI vs (Ic - Is)
                backgroundInt = patchNode.epochList.firstValue.protocolSettings('backgroundIntensity');
                discInt = patchNode.epochList.firstValue.protocolSettings('equivalentIntensity');
                surroundInt = backgroundInt + backgroundInt * parameterizedSurroundContrasts;
                IntDiff = discInt - surroundInt;
                IntDiff = IntDiff ./ backgroundInt; %normalize by image background
                
                tempDiff = (imageResponseMatrix.image.mean(ll,:) - imageResponseMatrix.disc.mean(ll,:));
                tempNLI = tempDiff ./ ...
                    (abs(imageResponseMatrix.image.mean(ll,:)) + abs(imageResponseMatrix.disc.mean(ll,:)));
                
                patchMean(patchCt) = (discInt- backgroundInt) / backgroundInt;
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
        end %for image

    end %for cell

    %NLI vs (Ic - Is)
%     binAndPlotEquallyPopulatedBins(AllIntDiff(:),AllNLI(:), 9, fig13)
    binAndPlotEquallyPopulatedBins(AllIntDiff(:),AllRespDiff(:), 9, fig13)
    

    disp('-------------------')
    disp([num2str(length(NLIvalues.naturalSurround)), ' patches'])
    disp([num2str(pp), ' cells'])

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

    if ~isempty(figureID)
        makeAxisStruct(fig2,['LEDmodSeg_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig3,['LEDmodS_NLIvsS_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig4,['LEDmodS_NLInat_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig6,['LEDmodS_DvsS_',figureID] ,'RFSurroundFigs')
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
        end
    end
end