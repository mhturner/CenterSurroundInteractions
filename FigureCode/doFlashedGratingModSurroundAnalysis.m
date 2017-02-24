function doFlashedGratingModSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    addParameter(ip,'exportFigs',true,@islogical);
    
    egTraceSurroundContrast = 0.75;
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    metric = ip.Results.metric;
    exportFigs = ip.Results.exportFigs;
    
    %eg Grating vs null
    figure; clf;
    fig2=gca;
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Grating center')
    set(get(fig2,'YLabel'),'String','No center')
    set(gcf, 'WindowStyle', 'docked')
    
    %eg NLI vs surround contrast, lines for different grating contrasts
    figure; clf;
    fig3=gca;
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Surround contrast')
    set(get(fig3,'YLabel'),'String','NLI')
    set(gcf, 'WindowStyle', 'docked')
    
    %eg cell traces, no surround
    figure; clf;
    fig4=gca;
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','Time (s)')
    set(get(fig4,'YLabel'),'String','Resp')
    set(gcf, 'WindowStyle', 'docked')
    
    %eg cell traces, with surround
    figure; clf;
    fig5=gca;
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','Time (s)')
    set(get(fig5,'YLabel'),'String','Resp')
    set(gcf, 'WindowStyle', 'docked')
    
    %Pop. NLI vs surround contrast
    figure; clf;
    fig6=gca;
    set(fig6,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig6,'XLabel'),'String','Surround contrast')
    set(get(fig6,'YLabel'),'String','Mean NLI')
    set(gcf, 'WindowStyle', 'docked')
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(gratingContrast)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    popNLIMatrix = [];
    for pp = 1:length(populationNodes)
        cellNode = populationNodes{pp};
        thisCellIsExample = 0;
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        
        noGratingContrasts = cellNode.children.length;
        noSurrounds = cellNode.children(1).children.length;
        
        clear gratingResponseMatrix
        gratingResponseMatrix.image.mean = nan(noGratingContrasts,noSurrounds);
        gratingResponseMatrix.image.err = nan(noGratingContrasts,noSurrounds);
        gratingResponseMatrix.disc.mean = nan(noGratingContrasts,noSurrounds);
        gratingResponseMatrix.disc.err = nan(noGratingContrasts,noSurrounds);
        NLImatrix = nan(noGratingContrasts,noSurrounds);
        surroundContrastValues = nan(1,noSurrounds);
        gratingContrastValues = nan(noGratingContrasts,1);
        for ii = 1:noGratingContrasts %for each gratingContrast
            gratingNode = cellNode.children(ii);
            gratingContrastValues(ii) = gratingNode.splitValue;
            for ss = 1:noSurrounds %for each surround
                surroundNode = gratingNode.children(ss);
                surroundContrastValues(ss) = surroundNode.splitValue;

                imageResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('image').epochList,recType);
                discResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('intensity').epochList,recType);

                gratingResponseMatrix.image.mean(ii,ss) = imageResp.(metric).mean;
                gratingResponseMatrix.image.err(ii,ss) = imageResp.(metric).SEM;

                gratingResponseMatrix.disc.mean(ii,ss) = discResp.(metric).mean;
                gratingResponseMatrix.disc.err(ii,ss) = discResp.(metric).SEM;

                newNLIvalue = (imageResp.(metric).mean - discResp.(metric).mean) ./ ...
                   (abs(imageResp.(metric).mean) + abs(discResp.(metric).mean));
                NLImatrix(ii,ss) = newNLIvalue;
                if gratingNode.custom.get('isExample')
                    
                    if (surroundContrastValues(ss) == egTraceSurroundContrast)
                        %add example traces, with target surround
                        imageTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,recType);
                        discTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,recType);
                        addLineToAxis(imageTrace.timeVector, imageTrace.mean,'egImage',fig5,'k','-','none')
                        addLineToAxis(discTrace.timeVector, discTrace.mean,'egDisc',fig5,[0.7 0.7 0.7],'-','none')
                    elseif (surroundContrastValues(ss) == 0)
                        %add example traces, no surround
                        imageTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,recType);
                        discTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,recType);
                        addLineToAxis(imageTrace.timeVector, imageTrace.mean,'egImage',fig4,'k','-','none')
                        addLineToAxis(discTrace.timeVector, discTrace.mean,'egDisc',fig4,[0.7 0.7 0.7],'-','none')
                    end
                end
                
            end %for surround

            if gratingNode.custom.get('isExample') %add to example figs
                thisCellIsExample = 1;
                % resp scatter plot:
                colors = pmkmp(length(surroundContrastValues));
                for ss = 1:length(surroundContrastValues)
                    currentSurround = surroundContrastValues(ss);
                    if currentSurround == 0
                        plotMarker = 'x';
                    else 
                        plotMarker = 'o';
                    end
                    
                    
                    addLineToAxis(gratingResponseMatrix.image.mean(ii,ss),...
                        gratingResponseMatrix.disc.mean(ii,ss),...
                        ['mean',num2str(ss)],fig2,colors(ss,:),'none',plotMarker)
                    tempX = gratingResponseMatrix.image.mean(ii,ss) + ...
                        [-gratingResponseMatrix.image.err(ii,ss), gratingResponseMatrix.image.err(ii,ss)];
                    tempY = [gratingResponseMatrix.disc.mean(ii,ss), gratingResponseMatrix.disc.mean(ii,ss)];
                    addLineToAxis(tempX, tempY,['errX',num2str(ss)],fig2,colors(ss,:),'-','none')
                    tempX = [gratingResponseMatrix.image.mean(ii,ss), gratingResponseMatrix.image.mean(ii,ss)];
                    tempY = gratingResponseMatrix.disc.mean(ii,ss) + ...
                        [-gratingResponseMatrix.disc.err(ii,ss), gratingResponseMatrix.disc.err(ii,ss)];
                    addLineToAxis(tempX, tempY,['errY',num2str(ss)],fig2,colors(ss,:),'-','none')
                end
                limUp = max([gratingResponseMatrix.image.mean(ii,:)', gratingResponseMatrix.disc.mean(ii,:)']);
                limDown = min([gratingResponseMatrix.image.mean(ii,:)', gratingResponseMatrix.disc.mean(ii,:)']);
                addLineToAxis([limDown limUp],[limDown limUp],...
                        'unity',fig2,'k','--','none')
            end %if example plotting
        end %for grate contrast
        if (thisCellIsExample)
            % NLI vs surround contrast plot:
            colors = repmat(linspace(0,0.8,length(gratingContrastValues))',[1 3]);
            for ii = 1:length(gratingContrastValues)
                addLineToAxis(surroundContrastValues, NLImatrix(ii,:),['NLI',num2str(ii)],fig3,colors(ii,:),'-','o')
            end
        end
        popNLIMatrix = cat(3,popNLIMatrix, NLImatrix);
    end %for cell
    %population NLI stats:
    meanPop = mean(popNLIMatrix,3);
    stdPop = std(popNLIMatrix,[],3);
    semPop = stdPop ./ sqrt(pp);
    
    colors = repmat(linspace(0,0.8,length(gratingContrastValues))',[1 3]);
    for ii = 1:length(gratingContrastValues)
        addLineToAxis(surroundContrastValues, meanPop(ii,:),['mean',num2str(ii)],fig6,colors(ii,:),'-','o')
        addLineToAxis(surroundContrastValues, meanPop(ii,:) - semPop(ii,:),['errDown',num2str(ii)],fig6,colors(ii,:),'--','none')
        addLineToAxis(surroundContrastValues, meanPop(ii,:) + semPop(ii,:),['errUp',num2str(ii)],fig6,colors(ii,:),'--','none')
    end

    if (exportFigs)
        figID = ['FGmod_',recType];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_NLI_',recType];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_TraceNoS_',recType];
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_TraceWithS_',recType];
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')

        figID = ['FGmod_pop_NLI_',recType];
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')
    end
end