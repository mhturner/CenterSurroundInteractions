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
    
    %Image vs Disc
    figure; clf;
    fig2=gca;
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Image center')
    set(get(fig2,'YLabel'),'String','Disc center')
    
    %NLI vs surround contrast
    figure; clf;
    fig3=gca;
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Surround contrast')
    set(get(fig3,'YLabel'),'String','Mean NLI')
    
    %NLI for no surround vs natural surround intensity
    figure; clf;
    fig4=gca;
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','NLI no surround')
    set(get(fig4,'YLabel'),'String','NLI with natural surround')
    
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
            imageResponseMatrix.image.mean = nan(noPatches,noSurrounds);
            imageResponseMatrix.image.err = nan(noPatches,noSurrounds);
            imageResponseMatrix.disc.mean = nan(noPatches,noSurrounds);
            imageResponseMatrix.disc.err = nan(noPatches,noSurrounds);
            surroundContrastValues = nan(noPatches,noSurrounds);
            for ll = 1:noPatches
               patchCt = patchCt + 1;
               patchNode = imageNode.children(ll);
               parameterizedSurroundContrasts = convertJavaArrayList(...
                        patchNode.epochList.firstValue.protocolSettings('surroundContrast'));
               for ss = 1:noSurrounds
                   currentSurroundContrast = patchNode.children(ss).splitValue;
                   
                   imageResp = getResponseAmplitudeStats(patchNode.children(ss).childBySplitValue('image').epochList,recType);
                   discResp = getResponseAmplitudeStats(patchNode.children(ss).childBySplitValue('intensity').epochList,recType);
                   
                   imageResponseMatrix.image.mean(ll,ss) = imageResp.(metric).mean;
                   imageResponseMatrix.image.err(ll,ss) = imageResp.(metric).SEM;
                   
                   imageResponseMatrix.disc.mean(ll,ss) = discResp.(metric).mean;
                   imageResponseMatrix.disc.err(ll,ss) = discResp.(metric).SEM;
                   surroundContrastValues(ll,ss) = currentSurroundContrast;
                   
                   newNLIvalue = (imageResp.(metric).mean - discResp.(metric).mean) ./ ...
                           (abs(imageResp.(metric).mean) + abs(discResp.(metric).mean));
                   if ismember(currentSurroundContrast,parameterizedSurroundContrasts)
                       NLIvalues.paramSurrounds(patchCt,ss) = newNLIvalue;
                   else
                       NLIvalues.naturalSurround(patchCt) = newNLIvalue;
                   end
               end %for surround
               
               if patchNode.custom.get('isExample') %add to example figs
                    patchLocation = str2num(patchNode.splitValue);
                    imageRes = getNaturalImagePatchFromLocation(patchLocation,imageNode.splitValue,'imageSize',[250 250]);
                    figure(30); imagesc(imageRes.images{1}); colormap(gray); axis image; axis off;
                    presentedContrasts = surroundContrastValues(ll,:);
                    colors = pmkmp(length(presentedContrasts));

                    for cc = 1:length(imageResponseMatrix.image.mean(ll,:))
                        currentSurround = surroundContrastValues(ll,cc);
                        if ismember(currentSurround,parameterizedSurroundContrasts)
                            plotMarker = 'o';
                        else %natural surround intensity
                            plotMarker = 'x';
                        end
                        
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
                    
                    limUp = max([imageResponseMatrix.image.mean(ll,:)', imageResponseMatrix.disc.mean(ll,:)']);
                    limDown = min([imageResponseMatrix.image.mean(ll,:)', imageResponseMatrix.disc.mean(ll,:)']);
                    addLineToAxis([limDown limUp],[limDown limUp],...
                            'unity',fig2,'k','--','none')

                end %if example plotting

            end %for patch location
        
        end %for image
        
        
    end %for cell
    
    disp(size(NLIvalues.paramSurrounds,1));
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
    addLineToAxis([0 1],[0 1],'unity',fig4,'k','--','none')
    
    if ~isempty(figureID)
        makeAxisStruct(fig2,['LEDmodSeg_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig3,['LEDmodS_NLIvsS_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig4,['LEDmodS_NLInat_',figureID] ,'RFSurroundFigs')
    end
end