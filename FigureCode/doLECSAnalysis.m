function doLECSAnalysis(node,varargin)
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
    
    %Linearity of surround. Eg ON and OFF cell
    figure; clf;
    fig2=gca;
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Center + image surround')
    set(get(fig2,'YLabel'),'String','Center + linear surround')
    
    % Linearity of center modulated by linear surround. Eg ON and OFF cell
    figure; clf;
    fig3=gca;
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Center image + linear surround')
    set(get(fig3,'YLabel'),'String','Center disc + linear surround')
    
    % Linearity of center modulated by image surround. Eg ON and OFF cell
    figure; clf;
    fig4=gca;
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','Center image + image surround')
    set(get(fig4,'YLabel'),'String','Center disc + image surround')
    

    centerStims = {'Image','none','Image','Equiv','none','Equiv','Image','Equiv'};
    surroundStims = {'none','Image','Image','none','Equiv','Equiv','Equiv','Image'};
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(imageName)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    
    imCt = 0;
    for pp = 1:length(populationNodes)
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        for ss = 1:cellNode.children.length %for each image
            imCt = imCt + 1;
            imageNode = cellNode.children(ss);
            noPatches = imageNode.children(1).children(1).children.length;
            tempResponseMatrix = nan(noPatches,length(centerStims));
            for cc = 1:length(centerStims) %for each stim pair
                stimPairNode = imageNode.childBySplitValue(centerStims{cc}).childBySplitValue(surroundStims{cc});
                for ii = 1:stimPairNode.children.length %for each patch location
                    newResp = getResponseAmplitudeStats(stimPairNode.children(ii).epochList,recType);
                    tempResponseMatrix(ii,cc) = newResp.(metric).mean;
                end
            end
            %population stats:
            temp = (tempResponseMatrix(:,1) - tempResponseMatrix(:,4)) ./ ...
                (abs(tempResponseMatrix(:,1)) + abs(tempResponseMatrix(:,4)));
            NLI.centerOnly(imCt) = nanmean(temp);
            
            temp = (tempResponseMatrix(:,7) - tempResponseMatrix(:,6)) ./ ...
                (abs(tempResponseMatrix(:,7)) + abs(tempResponseMatrix(:,6)));
            NLI.linearSurround(imCt) = nanmean(temp);
            
            temp = (tempResponseMatrix(:,3) - tempResponseMatrix(:,8)) ./ ...
                (abs(tempResponseMatrix(:,3)) + abs(tempResponseMatrix(:,8)));
            NLI.imageSurround(imCt) = nanmean(temp);

            if imageNode.custom.get('isExample') %add to example figs
                if strcmp(cellInfo.cellType,'ONparasol')
                    tempColor = [0.5 0.5 0.5];
                elseif strcmp(cellInfo.cellType,'OFFparasol')
                    tempColor = [0 0 0];
                end
                upLim = max(tempResponseMatrix(:));
                downLim = min(tempResponseMatrix(:));
                % surround linearity: Image center in both. Image vs
                %       annulus surround
                addLineToAxis(tempResponseMatrix(:,3),tempResponseMatrix(:,7),...
                    [cellInfo.cellType, '_', imageNode.splitValue],fig2,tempColor,'none','o')
                addLineToAxis([downLim upLim],[downLim upLim],...
                    'unity',fig2,'k','--','none')
                addLineToAxis(0,0,[cellInfo.cellType,'_',cellInfo.cellID],fig2,'k','none','none')
                
                colors = pmkmp(3);
                if strcmp(cellInfo.cellType,'OFFparasol')
                    % Linearity of center modulated by linear surround
                    addLineToAxis(tempResponseMatrix(:,1),tempResponseMatrix(:,4),...
                        [cellInfo.cellType, '_', imageNode.splitValue,'_center'],fig3,'k','none','o')
                    addLineToAxis(tempResponseMatrix(:,7),tempResponseMatrix(:,6),...
                        [cellInfo.cellType, '_', imageNode.splitValue,'_linSurround'],fig3,colors(2,:),'none','o')
                    addLineToAxis([downLim upLim],[downLim upLim],...
                        'unity',fig3,'k','--','none')
                    addLineToAxis(0,0,[cellInfo.cellType,'_',cellInfo.cellID],fig3,'k','none','none')

                    % Linearity of center modulated by image surround
                    addLineToAxis(tempResponseMatrix(:,1),tempResponseMatrix(:,4),...
                        [cellInfo.cellType, '_', imageNode.splitValue,'_center'],fig4,'k','none','o')
                    addLineToAxis(tempResponseMatrix(:,3),tempResponseMatrix(:,8),...
                        [cellInfo.cellType, '_', imageNode.splitValue,'_imSurround'],fig4,colors(3,:),'none','o')
                    addLineToAxis([downLim upLim],[downLim upLim],...
                        'unity',fig4,'k','--','none')
                    addLineToAxis(0,0,[cellInfo.cellType,'_',cellInfo.cellID],fig4,'k','none','none')
                end
            end

        end %for image
    end %for cell
    
    if ~isempty(figureID)
        makeAxisStruct(fig2,['LECS_linXY_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig3,['LECS_modLinearXY_',figureID] ,'RFSurroundFigs')
        makeAxisStruct(fig4,['LECS_modImageXY_',figureID] ,'RFSurroundFigs')
    end
end