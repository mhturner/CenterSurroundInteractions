function doFlashedGratingModSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    addParameter(ip,'exportFigs',true,@islogical);
    
    %+/- 0.75 for Off parasols
%     egTraceSurroundContrast1 = 0.75;
%     egTraceSurroundContrast2 = -0.75;
    egTraceSurroundContrast1 = -0.5;
    egTraceSurroundContrast2 = 0.5;
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    metric = ip.Results.metric;
    exportFigs = ip.Results.exportFigs;
    
    
    figure; clf; fig2=gca; initFig(fig2,'Grating center','No center') %eg Grating vs null
    figure; clf; fig3=gca; initFig(fig3,'Surround contrast','NLI') %eg NLI vs surround contrast, lines for different grating contrasts
    figure; clf; fig4=gca; initFig(fig4,'Time (s)','Resp') %eg cell traces, no surround
    figure; clf; fig5=gca; initFig(fig5,'Time (s)','Resp') %eg cell traces, with surround 1
    figure; clf; fig13=gca; initFig(fig13,'Time (s)','Resp') %eg cell traces, with surround 2
    figure; clf; fig6=gca; initFig(fig6,'Surround contrast','Mean NLI') %Pop. NLI vs surround contrast
    
    figure; clf; fig7=gca; initFig(fig7,'Time (s)','Trial') %Raster - s = 0, grating
    figure; clf; fig8=gca; initFig(fig8,'Time (s)','Trial') %Raster - s = 0, no grating
    figure; clf; fig9=gca; initFig(fig9,'Time (s)','Trial') %Raster - s = eg surround 1, grating
    figure; clf; fig10=gca; initFig(fig10,'Time (s)','Trial') %Raster - s = eg surround 1, no grating
    figure; clf; fig11=gca; initFig(fig11,'Time (s)','Trial') %Raster - s = eg surround 2, no grating
    figure; clf; fig12=gca; initFig(fig12,'Time (s)','Trial') %Raster - s = eg surround 2, no grating

    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(gratingContrast)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    popDiffMatrix = [];
    FlashedGratingModulatedSurround = [];
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
        DiffMatrix = nan(noGratingContrasts,noSurrounds);
        surroundContrastValues = nan(1,noSurrounds);
        gratingContrastValues = nan(noGratingContrasts,1);
        
        FlashedGratingModulatedSurround{pp} = cell(noGratingContrasts,noSurrounds);
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

                newDiff = (imageResp.(metric).mean - discResp.(metric).mean);
                
                % Pull out data to save:
                amp = surroundNode.epochList.firstValue.protocolSettings('amp');
                sampleRate = surroundNode.epochList.firstValue.protocolSettings('sampleRate'); %Hz
                baselineTime = surroundNode.epochList.firstValue.protocolSettings('preTime'); %msec
                baselinePoints = (baselineTime / 1e3) * sampleRate; %msec -> datapoints
                %grating data matrix:
                gratingResponse = riekesuite.getResponseMatrix(surroundNode.childBySplitValue('image').epochList,amp);
                baselines = mean(gratingResponse(:,1:baselinePoints),2); %baseline for each trial
                FlashedGratingModulatedSurround{pp}{ii,ss}.gratingResponse = gratingResponse - repmat(baselines,1,size(gratingResponse,2));
                if strcmp(recType,'extracellular')
                    %get spike binary matrices, too
                    tempTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,recType,'attachSpikeBinary',true);
                    FlashedGratingModulatedSurround{pp}{ii,ss}.gratingResponse_binary = tempTrace.binary;
                end
                
                %no grating data matrix:
                noGratingResponse = riekesuite.getResponseMatrix(surroundNode.childBySplitValue('intensity').epochList,amp);
                baselines = mean(noGratingResponse(:,1:baselinePoints),2); %baseline for each trial
                FlashedGratingModulatedSurround{pp}{ii,ss}.noGratingResponse = noGratingResponse - repmat(baselines,1,size(noGratingResponse,2));
                if strcmp(recType,'extracellular')
                    %get spike binary matrices, too
                    tempTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true);
                    FlashedGratingModulatedSurround{pp}{ii,ss}.noGratingResponse_binary = tempTrace.binary;
                end
                
                
                DiffMatrix(ii,ss) = newDiff;
                if gratingNode.custom.get('isExample')
                    addLineToAxis(0,0,cellInfo.cellID,fig2,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig3,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig4,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig5,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig7,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig8,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig9,'k','none','none')
                    addLineToAxis(0,0,cellInfo.cellID,fig10,'k','none','none')
                    if (surroundContrastValues(ss) == egTraceSurroundContrast1)
                        %add example traces, with target surround # 1
                        imageTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        discTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        addLineToAxis(imageTrace.timeVector, imageTrace.mean,'egImage',fig5,'k','-','none')
                        addLineToAxis(discTrace.timeVector, discTrace.mean,'egDisc',fig5,[0.7 0.7 0.7],'-','none')
                        
                        if strcmp(recType,'extracellular') %raster plot - eg surround 1
                            addRastersToFigure(imageTrace.binary,fig9)
                            addRastersToFigure(discTrace.binary,fig10)
                        end
                    elseif (surroundContrastValues(ss) == egTraceSurroundContrast2)
                        %add example traces, with target surround # 2
                        imageTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        discTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        addLineToAxis(imageTrace.timeVector, imageTrace.mean,'egImage',fig13,'k','-','none')
                        addLineToAxis(discTrace.timeVector, discTrace.mean,'egDisc',fig13,[0.7 0.7 0.7],'-','none')
                        
                        if strcmp(recType,'extracellular') %raster plot - eg surround 2
                            addRastersToFigure(imageTrace.binary,fig11)
                            addRastersToFigure(discTrace.binary,fig12)
                        end

                    elseif (surroundContrastValues(ss) == 0)
                        %add example traces, no surround
                        imageTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        discTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,recType,'attachSpikeBinary',true,'PSTHsigma',10);
                        addLineToAxis(imageTrace.timeVector, imageTrace.mean,'egImage',fig4,'k','-','none')
                        addLineToAxis(discTrace.timeVector, discTrace.mean,'egDisc',fig4,[0.7 0.7 0.7],'-','none')
                        
                        if strcmp(recType,'extracellular') %raster plot - no surround
                            addRastersToFigure(imageTrace.binary,fig7)
                            addRastersToFigure(discTrace.binary,fig8)
                        end %end if raster
                    end %end if eg surrounds
                end %end if example
                
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
                addLineToAxis(surroundContrastValues, DiffMatrix(ii,:),['NLI',num2str(ii)],fig3,colors(ii,:),'-','o')
            end
        end
        if size(DiffMatrix,1) < 4 %didn't get thru all 4 grating contrasts
            tempMat = DiffMatrix;
            DiffMatrix = [];
            targetCC = [0.25, 0.5, 0.75, 0.9];
            for cc = 1:4
                tempInd = find(targetCC(cc) == gratingContrastValues);
                if isempty(tempInd) 
                    DiffMatrix(cc,:) = nan(1,9);
                else
                    DiffMatrix(cc,:) = tempMat(tempInd,:);
                end
            end
        end
        GratingContrast{pp} = gratingContrastValues;
        SurroundContrast{pp} = surroundContrastValues;
        popDiffMatrix = cat(3,popDiffMatrix, DiffMatrix);
    end %for cell
    save(['FlashedGratingModulatedSurround_', cellInfo.cellType,'_', recType,'.mat'],'FlashedGratingModulatedSurround','GratingContrast','SurroundContrast')
    %population NLI stats:
    meanPop = nanmean(popDiffMatrix,3);
    stdPop = nanstd(popDiffMatrix,[],3);
    nPop = sum(~isnan(popDiffMatrix),3);
    
    semPop = stdPop ./ sqrt(nPop);
    
    gratingContrastValues = [0.25 0.5 0.75 0.9];
    colors = repmat(linspace(0,0.8,length(gratingContrastValues))',[1 3]);
    for ii = 1:length(gratingContrastValues)
        addLineToAxis(surroundContrastValues, meanPop(ii,:),['mean',num2str(ii)],fig6,colors(ii,:),'-','o')
        addLineToAxis(surroundContrastValues, meanPop(ii,:) - semPop(ii,:),['errDown',num2str(ii)],fig6,colors(ii,:),'--','none')
        addLineToAxis(surroundContrastValues, meanPop(ii,:) + semPop(ii,:),['errUp',num2str(ii)],fig6,colors(ii,:),'--','none')
    end

    if strcmp(cellInfo.cellType,'ONparasol')
        cellTypeTag = 'on';
    elseif strcmp(cellInfo.cellType,'OFFparasol')
        cellTypeTag = '';
    end
    
    if (exportFigs)
        figID = ['FGmod_',recType,cellTypeTag];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_NLI_',recType,cellTypeTag];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_TraceNoS_',recType,cellTypeTag];
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_TraceWithS_',recType,cellTypeTag];
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')

        figID = ['FGmod_pop_NLI_',recType,cellTypeTag];
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')
        
        figID = ['FGmod_TraceWithS2_',recType,cellTypeTag];
        makeAxisStruct(fig13,figID ,'RFSurroundFigs')
        
        if strcmp(recType,'extracellular') %raster plots
            figID = ['FGmod_rast_imAlone',cellTypeTag];
            makeAxisStruct(fig7,figID ,'RFSurroundFigs')
            
            figID = ['FGmod_rast_discAlone',cellTypeTag];
            makeAxisStruct(fig8,figID ,'RFSurroundFigs')
            
            figID = ['FGmod_rast_imSurr',cellTypeTag];
            makeAxisStruct(fig9,figID ,'RFSurroundFigs')
            
            figID = ['FGmod_rast_discSurr',cellTypeTag];
            makeAxisStruct(fig10,figID ,'RFSurroundFigs')
            
            figID = ['FGmod_rast_imSurr2',cellTypeTag];
            makeAxisStruct(fig11,figID ,'RFSurroundFigs')
            
            figID = ['FGmod_rast_discSurr2',cellTypeTag];
            makeAxisStruct(fig12,figID ,'RFSurroundFigs')
        
        end
    end
end