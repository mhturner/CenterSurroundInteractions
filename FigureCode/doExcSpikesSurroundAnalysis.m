function doExcSpikesSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    metric = ip.Results.metric;
    exportFigs = ip.Results.exportFigs;
    
        
    figure; clf; fig1=gca; initFig(fig1,'Time (s)','Response (spk/sec)') %spikes, traces
    figure; clf; fig2=gca; initFig(fig2,'Spot Size','Response (spikes)') %spikes, Area Summation
    
    figure; clf; fig3=gca; initFig(fig3,'Time (s)','Response (pA)') %exc, traces
    figure; clf; fig4=gca; initFig(fig4,'Spot Size','Response (pC)') %exc, Area Summation
    
    figure; clf; fig5=gca; initFig(fig5,'Center size, exc','Center size, spikes') % center size exc vs. spikes
    figure; clf; fig6=gca; initFig(fig6,'Surround strength, exc','Surround strength, spikes') %exc, Area Summation

    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                '@(list)splitOnRecKeyword(list)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    
    centerSize = [];
    surroundStrength = [];
    for pp = 1:length(populationNodes)
        currentNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(currentNode.epochList);
        % % % % % % % % % SPIKES % % % % % % % % % % % % % % % % % % % % % % % % 
        spikesNode = currentNode.childBySplitValue('extracellular');
        
        respAmps = nan(1,spikesNode.children.length);
        respErr = nan(1,spikesNode.children.length);
        spotSizes = nan(1,spikesNode.children.length);
        colors = pmkmp(spikesNode.children.length);
        for ee = 1:spikesNode.children.length %for spot sizes
            recType = getRecordingTypeFromEpochList(spikesNode.epochList);
            stats = getResponseAmplitudeStats(spikesNode.children(ee).epochList,recType);
            respAmps(ee) = stats.(metric).mean;
            respErr(ee) = stats.(metric).SEM;
            spotSizes(ee) = spikesNode.children(ee).splitValue;
            
            if currentNode.custom.get('isExample')
                responseTrace = getMeanResponseTrace(spikesNode.children(ee).epochList,recType,'PSTHsigma',10);
                addLineToAxis(responseTrace.timeVector,responseTrace.mean,...
                    ['spot',num2str(spotSizes(ee))],fig1,colors(ee,:),'-','none')
            end
        end %end for spot sizes
        
        %fit linear RF area-summation models:
        params0 = [max(respAmps), 40, max(respAmps), 150];
        [Kc,sigmaC,Ks,sigmaS] = fitDoGAreaSummation(spotSizes,respAmps,params0);
        fitX = 0:max(spotSizes);
        fitY = DoGAreaSummation([Kc,sigmaC,Ks,sigmaS], fitX);
        [~, ind] = max(fitY); 
        centerSize(pp,1) = fitX(ind);
        surroundStrength(pp,1) = (max(respAmps) - respAmps(end)) / max(respAmps);

        if currentNode.custom.get('isExample')
            addLineToAxis(spotSizes,respAmps,'data',fig2,'k','none','o')
            addLineToAxis(fitX,fitY,'fit',fig2,'k','-','none')
            addLineToAxis(spotSizes,respAmps + respErr,'errUp',fig2,'k','--','none')
            addLineToAxis(spotSizes,respAmps - respErr,'errDown',fig2,'k','--','none')
            
            addLineToAxis(0,0,cellInfo.cellID,fig1,'k','none','none')
            addLineToAxis(0,0,cellInfo.cellID,fig2,'k','none','none')
        end
        
        % % % % % % % % % EXC % % % % % % % % % % % % % % % % % % % % % % % % 
        excNode = currentNode.childBySplitValue('exc');
        respAmps = nan(1,excNode.children.length);
        respErr = nan(1,spikesNode.children.length);
        spotSizes = nan(1,excNode.children.length);
        colors = pmkmp(excNode.children.length);
        for ee = 1:excNode.children.length %for spot sizes
            recType = getRecordingTypeFromEpochList(excNode.epochList);
            stats = getResponseAmplitudeStats(excNode.children(ee).epochList,recType);
            respAmps(ee) = stats.(metric).mean;
            respErr(ee) = stats.(metric).SEM;
            spotSizes(ee) = excNode.children(ee).splitValue;
            
            if currentNode.custom.get('isExample')
                responseTrace = getMeanResponseTrace(excNode.children(ee).epochList,recType);
                addLineToAxis(responseTrace.timeVector,responseTrace.mean,...
                    ['spot',num2str(spotSizes(ee))],fig3,colors(ee,:),'-','none')
            end
        end %end for spot sizes
        
        %fit linear RF area-summation models:
        params0 = [max(respAmps), 40, max(respAmps), 150];
        [Kc,sigmaC,Ks,sigmaS] = fitDoGAreaSummation(spotSizes,respAmps,params0);
        fitX = 0:max(spotSizes);
        fitY = DoGAreaSummation([Kc,sigmaC,Ks,sigmaS], fitX);
        [~, ind] = max(fitY); 
        centerSize(pp,2) = fitX(ind);
        surroundStrength(pp,2) = (max(respAmps) - respAmps(end)) / max(respAmps);
        
        if currentNode.custom.get('isExample')
            addLineToAxis(spotSizes,respAmps,'data',fig4,'k','none','o')
            addLineToAxis(fitX,fitY,'fit',fig4,'k','-','none')
            addLineToAxis(spotSizes,respAmps + respErr,'errUp',fig4,'k','--','none')
            addLineToAxis(spotSizes,respAmps - respErr,'errDown',fig4,'k','--','none')
            
            addLineToAxis(0,0,cellInfo.cellID,fig3,'k','none','none')
            addLineToAxis(0,0,cellInfo.cellID,fig4,'k','none','none')
        end
        
    end %end for population
    addLineToAxis(centerSize(:,2),centerSize(:,1),'data',fig5,'k','none','o')
    addLineToAxis([0 max(centerSize(:))],[0 max(centerSize(:))],'unity',fig5,'k','--','none')
    
    addLineToAxis(surroundStrength(:,2),surroundStrength(:,1),'data',fig6,'k','none','o')
    addLineToAxis([0 max(surroundStrength(:))],[0 max(surroundStrength(:))],'unity',fig6,'k','--','none')
     
    if (exportFigs)
        figID = 'ES_spike_trace';
        makeAxisStruct(fig1,figID ,'RFSurroundFigs')
        
        figID = 'ES_spike_AS';
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')
        
        figID = 'ES_exc_trace';
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')
        
        figID = 'ES_exc_AS';
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')
        
        figID = 'ES_centerSize';
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')
        
        figID = 'ES_surroundStrength';
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')
    end
    
end