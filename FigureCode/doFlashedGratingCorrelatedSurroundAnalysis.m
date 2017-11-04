function doFlashedGratingCorrelatedSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    addParameter(ip,'exportFigs',true,@islogical);

    ip.parse(node,varargin{:});
    node = ip.Results.node;
    metric = ip.Results.metric;
    exportFigs = ip.Results.exportFigs;
    
    
    targetContrast = 0.5;
    targetIntensityValues = [0.3 0.4 0.5 0.6];
    
    egIntensityValue = 0.3;
        
    figure; clf; fig2=gca; initFig(fig2,'Center Intensity','NLI') % mean NLI vs. int for diff surrounds
    
    figure; clf; fig3=gca; initFig(fig3,'Time','trial') % Raster: grating, none
    figure; clf; fig4=gca; initFig(fig4,'Time','trial') % Raster: grating, corr
    figure; clf; fig5=gca; initFig(fig5,'Time','trial') % Raster: grating, acorr
 
    figure; clf; fig6=gca; initFig(fig6,'Time','trial') % Raster: disc, none
    figure; clf; fig7=gca; initFig(fig7,'Time','trial') % Raster: disc, corr
    figure; clf; fig8=gca; initFig(fig8,'Time','trial') % Raster: disc, acorr
    
    figure; clf; fig9=gca; initFig(fig9,'no mean, add mean, add surround','resp. diff. (spk)') % resp. diff by category, lines
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(gratingContrast)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    nliMat = [];
    diffMat = [];

    for pp = 1:length(populationNodes)
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        
        gratingNode = cellNode.childBySplitValue(targetContrast);
            
        respMat = nan(length(targetIntensityValues), 3, 2); %intensity x surroundType x image/disc;
        for ii = 1:length(targetIntensityValues)%for each currentIntensity
            
            intensityNode = gratingNode.childBySplitValue(targetIntensityValues(ii));

            for ss = 1:3 %for each surround type (acorr, corr, none)
                switch ss
                    case 1
                        surroundNode = intensityNode.childBySplitValue('acorr');
                    case 2
                        surroundNode = intensityNode.childBySplitValue('corr');
                    case 3
                        surroundNode = intensityNode.childBySplitValue('none');
                end
                

                imageResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('image').epochList,recType);
                discResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('intensity').epochList,recType);
                
                if cellNode.custom.get('isExample') & (intensityNode.splitValue == egIntensityValue) %#ok<AND2>
                    imageTrace = getMeanResponseTrace(surroundNode.childBySplitValue('image').epochList,...
                        recType,'attachSpikeBinary',true,'PSTHsigma',10);
                    discTrace = getMeanResponseTrace(surroundNode.childBySplitValue('intensity').epochList,...
                        recType,'attachSpikeBinary',true,'PSTHsigma',10);
                    switch ss
                        case 1 %acorr
                            addRastersToFigure(imageTrace.binary,fig5)
                            addRastersToFigure(discTrace.binary,fig8)
                        case 2 %corr
                            addRastersToFigure(imageTrace.binary,fig4)
                            addRastersToFigure(discTrace.binary,fig7)
                        case 3 %none
                            addRastersToFigure(imageTrace.binary,fig3)
                            addRastersToFigure(discTrace.binary,fig6)
                    end
                    
                end

                respMat(ii,ss,1) = imageResp.(metric).mean;
                respMat(ii,ss,2) = discResp.(metric).mean;

            end %end for surroundTag
        end %end for currentIntensity      
        
        diffMat(:,:,pp) = (respMat(:,:,1) - respMat(:,:,2)); %#ok<AGROW>
%         nliMat(:,:,pp) = (respMat(:,:,1) - respMat(:,:,2)) ./ (respMat(:,:,1) + respMat(:,:,2)); %#ok<AGROW>
tempMat = [respMat(:,:,1), respMat(:,:,2)];
normFactor = mean(tempMat(:));
        nliMat(:,:,pp) = (respMat(:,:,1) - respMat(:,:,2)) ./ normFactor; %#ok<AGROW>
    end %end for cell

    
    useMat = nliMat;
    intDrawInd = 1; %1 = 0.3, strongest negative mean
    %means, errs in categories
    noMean.meanResp = mean(useMat(3,3,:),3);
    noMean.errResp = std(useMat(3,3,:),[],3) ./ sqrt(size(useMat,3));
    
    withMean.meanResp = mean(useMat(intDrawInd,3,:),3);
    withMean.errResp = std(useMat(intDrawInd,3,:),[],3) ./ sqrt(size(useMat,3));
    
    corrS.meanResp = mean(useMat(intDrawInd,2,:),3);
    corrS.errResp = std(useMat(intDrawInd,2,:),[],3) ./ sqrt(size(useMat,3));
    
    acorrS.meanResp = mean(useMat(intDrawInd,1,:),3);
    acorrS.errResp = std(useMat(intDrawInd,1,:),[],3) ./ sqrt(size(useMat,3));
    
    addLineToAxis(1, noMean.meanResp,'noMean_mean',fig9,'k','none','s')
    addLineToAxis([1, 1], noMean.meanResp + [noMean.errResp, -noMean.errResp],...
        'noMean_err',fig9,'k','-','none')
    
    addLineToAxis(2, withMean.meanResp,'withMean_mean',fig9,'k','none','s')
    addLineToAxis([2, 2], withMean.meanResp + [withMean.errResp, -withMean.errResp],...
        'withMean_err',fig9,'k','-','none')
    
    addLineToAxis(3, corrS.meanResp,'corrS_mean',fig9,'g','none','s')
    addLineToAxis([3, 3], corrS.meanResp + [corrS.errResp, -corrS.errResp],...
        'corrS_err',fig9,'g','-','none')
    
    addLineToAxis(3, acorrS.meanResp,'acorrS_mean',fig9,'r','none','s')
    addLineToAxis([3, 3], acorrS.meanResp + [acorrS.errResp, -acorrS.errResp],...
        'acorrS_err',fig9,'r','-','none')
    
    [~, p] = ttest(useMat(3,3,:), useMat(intDrawInd,3,:)) %addition of mean
    [~, p] = ttest(useMat(intDrawInd,3,:), useMat(intDrawInd,2,:)) %addition of corr S
    [~, p] = ttest(useMat(intDrawInd,3,:), useMat(intDrawInd,1,:)) %addition of acorr S
    [~, p] = ttest(useMat(intDrawInd,2,:), useMat(intDrawInd,1,:)) %compare corr and acorr s
    for cc = 1:size(useMat,3)
        %lines
        addMeanLine = [useMat(3,3,cc), useMat(intDrawInd,3,cc)];
        addLineToAxis([1,2], addMeanLine,['addMeanLine',num2str(cc)],fig9,'k','-','o')

        corrLine = [useMat(intDrawInd,3,cc), useMat(intDrawInd,2,cc)];
        addLineToAxis([2,3], corrLine,['corrLine',num2str(cc)],fig9,'g','-','o')

        acorrLine = [useMat(intDrawInd,3,cc), useMat(intDrawInd,1,cc)];
        addLineToAxis([2,3], acorrLine,['acorrLine',num2str(cc)],fig9,'r','-','o')
    end
    
    %NLI stat calcs:
    nMat = sum(~isnan(nliMat),3);
    meanMat = nanmean(nliMat,3);
    errMat = nanstd(nliMat,[],3) ./ sqrt(nMat);

    addLineToAxis(targetIntensityValues,meanMat(:,1),'acorr',fig2,'r','-','o')
    addLineToAxis(targetIntensityValues,meanMat(:,1) + errMat(:,1),'acorr_eUp',fig2,'r','--','none')
    addLineToAxis(targetIntensityValues,meanMat(:,1) - errMat(:,1),'acorr_eDown',fig2,'r','--','none')

    addLineToAxis(targetIntensityValues,meanMat(:,2),'corr',fig2,'g','-','o')
    addLineToAxis(targetIntensityValues,meanMat(:,2) + errMat(:,2),'corr_eUp',fig2,'g','--','none')
    addLineToAxis(targetIntensityValues,meanMat(:,2) - errMat(:,2),'corr_eDown',fig2,'g','--','none')
    
    addLineToAxis(targetIntensityValues,meanMat(:,3),'none',fig2,'k','-','o')
    addLineToAxis(targetIntensityValues,meanMat(:,3) + errMat(:,3),'none_eUp',fig2,'k','--','none')
    addLineToAxis(targetIntensityValues,meanMat(:,3) - errMat(:,3),'none_eDown',fig2,'k','--','none')
    
    
    if (exportFigs)
        makeAxisStruct(fig2,'FGcorrS_summary' ,'RFSurroundFigs')
        
        makeAxisStruct(fig3,'FGcorrS_rast_gn' ,'RFSurroundFigs')
        makeAxisStruct(fig4,'FGcorrS_rast_gc' ,'RFSurroundFigs')
        makeAxisStruct(fig5,'FGcorrS_rast_ga' ,'RFSurroundFigs')
        
        makeAxisStruct(fig6,'FGcorrS_rast_dn' ,'RFSurroundFigs')
        makeAxisStruct(fig7,'FGcorrS_rast_dc' ,'RFSurroundFigs')
        makeAxisStruct(fig8,'FGcorrS_rast_da' ,'RFSurroundFigs')
        
        makeAxisStruct(fig9,'FGcorrS_category' ,'RFSurroundFigs')

    end
end