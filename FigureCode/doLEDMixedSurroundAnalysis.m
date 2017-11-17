function doLEDMixedSurroundAnalysis(node,varargin)
    ip = inputParser;
    expectedMetrics = {'integrated','peak'};
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'metric','integrated',...
        @(x) any(validatestring(x,expectedMetrics)));
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    metric = ip.Results.metric;
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(imageName)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    
    figure; clf; fig2=gca; initFig(fig2,'Diff, none (spk)','Diff, nat (spk)') % image mean resp diff, none Vs nat
    figure; clf; fig3=gca; initFig(fig3,'Diff, none (spk)','Diff, mix (spk)') % image mean resp diff, none Vs mix
    figure; clf; fig4=gca; initFig(fig4,'Diff, nat (spk)','Diff, mix (spk)') % image mean resp diff, none Vs mix
    
    figure; clf; fig5=gca; initFig(fig5,'Diff, none (spk)','Diff, nat (spk)') % cell mean resp diff, none Vs nat
    figure; clf; fig6=gca; initFig(fig6,'Diff, none (spk)','Diff, mix (spk)') % cell mean resp diff, none Vs mix
    figure; clf; fig7=gca; initFig(fig7,'Diff, nat (spk)','Diff, mix (spk)') % cell mean resp diff, none Vs mix
    
    figure; clf; fig8=gca; initFig(fig8,'Resp. to image (spk)','Resp. to disc (spk)') % eg cell: IvD - none
    figure; clf; fig9=gca; initFig(fig9,'Resp. to image (spk)','Resp. to disc (spk)') % eg cell: IvD - nat
    figure; clf; fig10=gca; initFig(fig10,'Resp. to image (spk)','Resp. to disc (spk)') % eg cell: IvD - mix
    

    allImageResponseMatrix = [];
    allDiscResponseMatrix = [];
    allImageIDs = [];
    allCellIDs = [];
    allEqContrast = [];
    allSpatialContrast = [];
    allMixSurroundInt = [];
    allBackgroundInt = [];

    imageCt = 0;
    cellCt = 0;
    
    filterSize = 106;
    [rr, cc] = meshgrid(1:filterSize,1:filterSize);
    surroundBinary = logical(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=103 & ...
        sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)>18);
    for pp = 1:length(populationNodes) % for cell in pop
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        cellCt = cellCt + 1;
        
        ImageResponseMatrix = []; % rows = center image, 3 columns = surround condition ([none, natural, mixed])
        DiscResponseMatrix = [];
        eqContrast = [];
        spatialContrast = [];
        backgroundInt = [];
        mixSurroundInt = [];
        for imageInd = 1:cellNode.children.length % for image ID
            imageCt = imageCt + 1;
            imageNode = cellNode.children(imageInd);
            for patchInd = 1:imageNode.children.length % for center patch
                patchNode = imageNode.children(patchInd);    
                eqInt = patchNode.epochList.firstValue.protocolSettings('equivalentIntensity');
                bg = patchNode.epochList.firstValue.protocolSettings('backgroundIntensity');
                eqContrast(patchInd) = (eqInt - bg) / bg;
                backgroundInt(patchInd) = bg;
                
                patchLoc = convertJavaArrayList(patchNode.epochList.firstValue.protocolSettings('currentCenterLocation'));
                imageName = patchNode.epochList.firstValue.protocolSettings('imageName');
                imageRes = getNaturalImagePatchFromLocation(patchLoc,imageName,'imageSize',[300 300]);
                %spatialContrast = stdev / mean
                spatialContrast(patchInd) = std(imageRes.images{1}(:)) / mean(imageRes.images{1}(:));
                mixLoc = convertJavaArrayList(patchNode.epochList.firstValue.protocolSettings('currentMixedSurroundPatchLocation'));
                mixRes = getNaturalImagePatchFromLocation(mixLoc,imageName,'imageSize',[700 700]);
                mixSurroundInt(patchInd) = mean(mixRes.images{1}(surroundBinary));

                %get image and disc responses, no surround:
                noSurroundNode = patchNode.childBySplitValue('none');
                imageResp = getResponseAmplitudeStats(noSurroundNode.childBySplitValue('image').epochList,recType);
                discResp = getResponseAmplitudeStats(noSurroundNode.childBySplitValue('intensity').epochList,recType);
                ImageResponseMatrix(patchInd,1) = imageResp.(metric).mean;
                DiscResponseMatrix(patchInd,1) = discResp.(metric).mean;
                
                
                %get image and disc responses, nat surround:
                natSurroundNode = patchNode.childBySplitValue('nat');
                imageResp = getResponseAmplitudeStats(natSurroundNode.childBySplitValue('image').epochList,recType);
                discResp = getResponseAmplitudeStats(natSurroundNode.childBySplitValue('intensity').epochList,recType);
                ImageResponseMatrix(patchInd,2) = imageResp.(metric).mean;
                DiscResponseMatrix(patchInd,2) = discResp.(metric).mean;
                
                
                %get image and disc responses, mixed surround:
                mixedSurroundNode = patchNode.childBySplitValue('mixed');
                imageResp = getResponseAmplitudeStats(mixedSurroundNode.childBySplitValue('image').epochList,recType);
                discResp = getResponseAmplitudeStats(mixedSurroundNode.childBySplitValue('intensity').epochList,recType);
                ImageResponseMatrix(patchInd,3) = imageResp.(metric).mean;
                DiscResponseMatrix(patchInd,3) = discResp.(metric).mean;
            end % for center patch
            
            if imageNode.custom.get('isExample')
                unityUp = 1.1*max([ImageResponseMatrix(:), DiscResponseMatrix(:)]);
                addLineToAxis(ImageResponseMatrix(:,1),DiscResponseMatrix(:,1),'data',fig8,'k','none','o')
                addLineToAxis([0 unityUp],[0 unityUp],'unity',fig8,'k','--','none')

                addLineToAxis(ImageResponseMatrix(:,2),DiscResponseMatrix(:,2),'data',fig9,'g','none','o')
                addLineToAxis([0 unityUp],[0 unityUp],'unity',fig9,'k','--','none')

                addLineToAxis(ImageResponseMatrix(:,3),DiscResponseMatrix(:,3),'data',fig10,'r','none','o')
                addLineToAxis([0 unityUp],[0 unityUp],'unity',fig10,'k','--','none')
            end
            
            %compute NLI matrix for this image
            % rows = center patch, 3 columns = surround condition ([none, natural, mixed])
            allImageResponseMatrix = cat(1,allImageResponseMatrix,ImageResponseMatrix);
            allDiscResponseMatrix = cat(1,allDiscResponseMatrix,DiscResponseMatrix);
            %cell/image index each patch belongs to, for within-cell/image stats:
            noPatches = size(ImageResponseMatrix,1); %no patches in this image
            allImageIDs = cat(2,allImageIDs,imageCt.*ones(1,noPatches));
            allCellIDs = cat(2,allCellIDs,cellCt.*ones(1,noPatches));
            %image patch stats:
            allEqContrast = cat(2,allEqContrast,eqContrast);
            allSpatialContrast = cat(2,allSpatialContrast,spatialContrast);
            allMixSurroundInt = cat(2,allMixSurroundInt,mixSurroundInt);
            allBackgroundInt = cat(2,allBackgroundInt,backgroundInt);
        end % for image ID
    end % for cell in pop

    allDiffMatrix = allImageResponseMatrix - allDiscResponseMatrix;
    allNliMatrix = allDiffMatrix ./ (allImageResponseMatrix + allDiscResponseMatrix);
% % % % % % MIX S: NLI VS Ic - Is % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    allCenterInt = (allEqContrast .* allBackgroundInt) + allBackgroundInt;
    relativeCSdiff = (allCenterInt - allMixSurroundInt) ./ allBackgroundInt;
    tempMixVex = allNliMatrix(:,3);

    noBins = 9;
    figure; clf; fig12=gca; initFig(fig12,'Ic - Is','NLI')
    binAndPlotEquallyPopulatedBins(relativeCSdiff(~isnan(tempMixVex)),tempMixVex(~isnan(tempMixVex)),noBins,fig12,[0.5 0.5 0.5],'mix')

% % % % % % NLI VS. BINNED CENTER MEAN % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    noBins = 6;
    cutInds = find(allEqContrast > 0); %cut positive means, too few spike responses
    allEqContrast(cutInds) = [];
    allNliMatrix(cutInds,:) = [];

    figure; clf; fig11=gca; initFig(fig11,'Relative center mean','NLI') % NLI vs center intensity
    tempNLIs = allNliMatrix(:,1);
    binAndPlotEquallyPopulatedBins(allEqContrast(~isnan(tempNLIs)),tempNLIs(~isnan(tempNLIs)),noBins,fig11,'k','none')

    tempNLIs = allNliMatrix(:,2);
    binAndPlotEquallyPopulatedBins(allEqContrast(~isnan(tempNLIs)),tempNLIs(~isnan(tempNLIs)),noBins,fig11,'g','nat')

    tempNLIs = allNliMatrix(:,3);
    binAndPlotEquallyPopulatedBins(allEqContrast(~isnan(tempNLIs)),tempNLIs(~isnan(tempNLIs)),noBins,fig11,'r','mix')
% % % % % % MEAN BY IMAGE % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    colors = hsv(cellCt);
    imageMeans = [];
    for ii = 1:imageCt
        pullInds = find(allImageIDs == ii);
        cellInd = allCellIDs(pullInds(1));
        imageMeans(ii,:) = [mean(allDiffMatrix(pullInds,1)),mean(allDiffMatrix(pullInds,2)),mean(allDiffMatrix(pullInds,3))];
        
        addLineToAxis(mean(allDiffMatrix(pullInds,1)),mean(allDiffMatrix(pullInds,2)),['i',num2str(ii)],fig2,colors(cellInd,:),'none','o')
        addLineToAxis(mean(allDiffMatrix(pullInds,1)),mean(allDiffMatrix(pullInds,3)),['i',num2str(ii)],fig3,colors(cellInd,:),'none','o')
        addLineToAxis(mean(allDiffMatrix(pullInds,2)),mean(allDiffMatrix(pullInds,3)),['i',num2str(ii)],fig4,colors(cellInd,:),'none','o')
    end
    %stats:
    %paired t-test
    [~, p] = ttest(imageMeans(:,1),imageMeans(:,2)) %none, nat
    [~, p] = ttest(imageMeans(:,1),imageMeans(:,3)) %none, mix
    [~, p] = ttest(imageMeans(:,2),imageMeans(:,3)) %nat, mix
    
    %mean & error points:
    tempMean = mean(imageMeans);
    tempErr = std(imageMeans) ./ sqrt(imageCt);
    
    addLineToAxis(tempMean(1),tempMean(2),'mean',fig2,'k','none','.')
    addLineToAxis(tempMean(1) + tempErr(1).*[-1,1],[tempMean(2),tempMean(2)],'errX',fig2,'k','-','none')
    addLineToAxis([tempMean(1), tempMean(1)],tempMean(2) + tempErr(2).*[-1,1],'errY',fig2,'k','-','none')
    
    addLineToAxis(tempMean(1),tempMean(3),'mean',fig3,'k','none','.')
    addLineToAxis(tempMean(1) + tempErr(1).*[-1,1],[tempMean(3),tempMean(3)],'errX',fig3,'k','-','none')
    addLineToAxis([tempMean(1), tempMean(1)],tempMean(3) + tempErr(3).*[-1,1],'errY',fig3,'k','-','none')
    
    addLineToAxis(tempMean(2),tempMean(3),'mean',fig4,'k','none','.')
    addLineToAxis(tempMean(2) + tempErr(2).*[-1,1],[tempMean(3),tempMean(3)],'errX',fig4,'k','-','none')
    addLineToAxis([tempMean(2), tempMean(2)],tempMean(3) + tempErr(3).*[-1,1],'errY',fig4,'k','-','none')
    
    addLineToAxis([0 5],[0 5],'unity',fig2,'k','--','none')
    addLineToAxis([0 5],[0 5],'unity',fig3,'k','--','none')
    addLineToAxis([0 5],[0 5],'unity',fig4,'k','--','none')
    
    
% % % % % % MEAN BY CELL % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    for ii = 1:cellCt
        pullInds = find(allCellIDs == ii);
        for pp = 1:3
            switch pp
                case 1
                    xInd = 1; yInd = 2; figID = fig5; color = 'g';
                case 2
                    xInd = 1; yInd = 3; figID = fig6; color = 'r';
                case 3
                    xInd = 2; yInd = 3; figID = fig7; color = 'k';
            end

            tempMeanX = mean(allDiffMatrix(pullInds,xInd)); tempMeanY = mean(allDiffMatrix(pullInds,yInd));
            tempErrX = std(allDiffMatrix(pullInds,xInd)) ./ sqrt(length(pullInds));
            tempErrY = std(allDiffMatrix(pullInds,yInd)) ./ sqrt(length(pullInds));
            addLineToAxis(tempMeanX,tempMeanY,['mean_',num2str(ii)],figID,color,'none','o')
            addLineToAxis([tempMeanX tempMeanX],tempMeanY + tempErrY.*[-1, 1],['yErr_',num2str(ii)],figID,color,'-','none')
            addLineToAxis(tempMeanX + tempErrX.*[-1, 1],[tempMeanY tempMeanY],['xErr_',num2str(ii)],figID,color,'-','none')
        end

    end
    addLineToAxis([0 5],[0 5],'unity',fig5,'k','--','none')
    addLineToAxis([0 5],[0 5],'unity',fig6,'k','--','none')
    addLineToAxis([0 5],[0 5],'unity',fig7,'k','--','none')

% % % % % % POPULATION N INFO % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    disp('No. patches:')
    disp(size(allImageResponseMatrix,1));
    disp('No. images:')
    disp(imageCt);
    disp('No. cells:')
    disp(cellCt);
    drawnow();

    if (exportFigs)
        makeAxisStruct(fig2,'MixS_nat' ,'RFSurroundFigs')
        makeAxisStruct(fig3,'MixS_mix' ,'RFSurroundFigs')
        makeAxisStruct(fig4,'MixS_natmix' ,'RFSurroundFigs')

        makeAxisStruct(fig5,'MixS_natC' ,'RFSurroundFigs')
        makeAxisStruct(fig6,'MixS_mixC' ,'RFSurroundFigs')
        makeAxisStruct(fig7,'MixS_natmixC' ,'RFSurroundFigs')

        makeAxisStruct(fig8,'MixS_egNone' ,'RFSurroundFigs')
        makeAxisStruct(fig9,'MixS_egNat' ,'RFSurroundFigs')
        makeAxisStruct(fig10,'MixS_egMix' ,'RFSurroundFigs')
        
        makeAxisStruct(fig11,'MixS_NLIvCtr' ,'RFSurroundFigs')
        makeAxisStruct(fig12,'MixS_NLIvMixCS' ,'RFSurroundFigs')
    end
    
end