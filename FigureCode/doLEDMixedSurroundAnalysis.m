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
    
    % Image vs. disc - no surround
    figure; clf;
    fig6=gca;
    set(fig6,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 14)
    set(get(fig6,'XLabel'),'String','Resp. to image (spk)')
    set(get(fig6,'YLabel'),'String','Resp. to disc (spk)')
    set(gcf, 'WindowStyle', 'docked')
    
    % Image vs. disc - natural surround
    figure; clf;
    fig7=gca;
    set(fig7,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 14)
    set(get(fig7,'XLabel'),'String','Resp. to image (spk)')
    set(get(fig7,'YLabel'),'String','Resp. to disc (spk)')
    set(gcf, 'WindowStyle', 'docked')
    
    % Image vs. disc - mixed surround
    figure; clf;
    fig8=gca;
    set(fig8,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 14)
    set(get(fig8,'XLabel'),'String','Resp. to image (spk)')
    set(get(fig8,'YLabel'),'String','Resp. to disc (spk)')
    set(gcf, 'WindowStyle', 'docked')
    
    % NLI vs mean center intensity
    figure; clf;
    fig9=gca;
    set(fig9,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 14)
    set(get(fig9,'XLabel'),'String','Center contrast')
    set(get(fig9,'YLabel'),'String','NLI')
    set(gcf, 'WindowStyle', 'docked')
    
    

    meanNLI = [];
    allNLIMatrix = [];
    
    meanDiff = [];
    allDiffMatrix = [];
    
    allEqContrast = [];
    allSpatialContrast = [];
    allImageResponseMatrix = [];
    
    p_none = [];
    p_nat = [];
    p_mix = [];
    imageCt = 0;
    for pp = 1:length(populationNodes) % for cell in pop
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        
        ImageResponseMatrix = []; % rows = center image, 3 columns = surround condition ([none, natural, mixed])
        DiscResponseMatrix = [];
        eqContrast = [];
        spatialContrast = [];
        for imageInd = 1:cellNode.children.length % for image ID
            imageCt = imageCt + 1;
            imageNode = cellNode.children(imageInd);
            for patchInd = 1:imageNode.children.length % for center patch
                patchNode = imageNode.children(patchInd);    
                eqInt = patchNode.epochList.firstValue.protocolSettings('equivalentIntensity');
                bg = patchNode.epochList.firstValue.protocolSettings('backgroundIntensity');
                eqContrast(patchInd) = (eqInt - bg) / bg;
                
                patchLoc = convertJavaArrayList(patchNode.epochList.firstValue.protocolSettings('currentCenterLocation'));
                imageName = patchNode.epochList.firstValue.protocolSettings('imageName');
                imageRes = getNaturalImagePatchFromLocation(patchLoc,imageName,'imageSize',[300 300]);
                %spatialContrast = stdev / mean
                spatialContrast(patchInd) = std(imageRes.images{1}(:)) / mean(imageRes.images{1}(:));

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
                addLineToAxis(ImageResponseMatrix(:,1),DiscResponseMatrix(:,1),'data',fig6,'k','none','o')
                addLineToAxis([0 unityUp],[0 unityUp],'unity',fig6,'k','--','none')

                addLineToAxis(ImageResponseMatrix(:,2),DiscResponseMatrix(:,2),'data',fig7,'g','none','o')
                addLineToAxis([0 unityUp],[0 unityUp],'unity',fig7,'k','--','none')

                addLineToAxis(ImageResponseMatrix(:,3),DiscResponseMatrix(:,3),'data',fig8,'r','none','o')
                addLineToAxis([0 unityUp],[0 unityUp],'unity',fig8,'k','--','none')
            end
            
            %compute NLI matrix for this image
            % rows = center patch, 3 columns = surround condition ([none, natural, mixed])
            diffMatrix = ImageResponseMatrix - DiscResponseMatrix;
            NLIresultsMatrix = diffMatrix ./ (abs(ImageResponseMatrix) + abs(DiscResponseMatrix));
            
            allDiffMatrix = cat(1,allDiffMatrix,diffMatrix);
            meanDiff(imageCt,:) = nanmean(diffMatrix,1);
            
            allNLIMatrix = cat(1,allNLIMatrix,NLIresultsMatrix);
            meanNLI(imageCt,:) = nanmean(NLIresultsMatrix,1);
            
            allImageResponseMatrix = cat(1,allImageResponseMatrix,ImageResponseMatrix);
            
            allEqContrast = cat(2,allEqContrast,eqContrast);
            allSpatialContrast = cat(2,allSpatialContrast,spatialContrast);
            
            
            %R2 values:
            %no surround
            imageResp = ImageResponseMatrix(:,1); discResp = DiscResponseMatrix(:,1);
            ssTot=sum((imageResp-mean(imageResp)).^2); 
            ssErr=sum((imageResp - discResp).^2);
            rSquared=1-ssErr/ssTot;
            p_none(imageCt) = rSquared;
            
            %nat surround
            imageResp = ImageResponseMatrix(:,2); discResp = DiscResponseMatrix(:,2);
            ssTot=sum((imageResp-mean(imageResp)).^2); 
            ssErr=sum((imageResp - discResp).^2);
            rSquared=1-ssErr/ssTot;
            p_nat(imageCt) = rSquared;
            
            %mix surround
            imageResp = ImageResponseMatrix(:,3:end); discResp = DiscResponseMatrix(:,3:end);
            imageResp = imageResp (:); discResp = discResp(:);
            ssTot=sum((imageResp-mean(imageResp)).^2); 
            ssErr=sum((imageResp - discResp).^2);
            rSquared=1-ssErr/ssTot;
            p_mix(imageCt) = rSquared;
            
            
        end % for image ID
    end % for cell in pop

    % response vs mean center intensity
    figure; clf;
    fig11=gca;
    set(fig11,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 14)
    set(get(fig11,'XLabel'),'String','Relative center intensity')
    set(get(fig11,'YLabel'),'String','Response (norm)')
    set(gcf, 'WindowStyle', 'docked')
    
    noBins = 12;
    [~,binMean,~,binID] = histcounts_equallyPopulatedBins(allEqContrast,noBins);
    
    highContrastInds = find(allSpatialContrast < median(allSpatialContrast));
    
    resMean = struct;
    resErr = struct;
    for bb = 1:noBins
        tempInds = find(binID == bb);
        pullInds = intersect(tempInds,highContrastInds);
        resMean.none(bb) = mean(allImageResponseMatrix(pullInds,1));
        resMean.nat(bb) = mean(allImageResponseMatrix(pullInds,2));
        resMean.mix(bb) = mean(allImageResponseMatrix(pullInds,3));
        
        resErr.none(bb) = std(allImageResponseMatrix(pullInds,1)) ./ sqrt(length(pullInds));
        resErr.nat(bb) = std(allImageResponseMatrix(pullInds,2)) ./ sqrt(length(pullInds));
        resErr.mix(bb) = std(allImageResponseMatrix(pullInds,3)) ./ sqrt(length(pullInds));
    end
    tempMean = resMean.none ./max(resMean.none);
    tempErr = resErr.none ./max(resMean.none);
    addLineToAxis(binMean,tempMean,'none',fig11,'k','-','.')
    addLineToAxis(binMean,tempMean + tempErr,'none_eUp',fig11,'k','--','none')
    addLineToAxis(binMean,tempMean - tempErr,'none_eDown',fig11,'k','--','none')
    
    tempMean = resMean.nat ./max(resMean.nat);
    tempErr = resErr.nat ./max(resMean.nat);
    addLineToAxis(binMean,tempMean,'nat',fig11,'g','-','.')
    addLineToAxis(binMean,tempMean + tempErr,'nat_eUp',fig11,'g','--','none')
    addLineToAxis(binMean,tempMean - tempErr,'nat_eDown',fig11,'g','--','none')
    
%     tempMean = resMean.mix ./max(resMean.mix);
%     tempErr = resErr.mix ./max(resMean.mix);
%     addLineToAxis(binMean,tempMean,'mix',fig11,'r','-','.')
%     addLineToAxis(binMean,tempMean + tempErr,'mix_eUp',fig11,'r','--','none')
%     addLineToAxis(binMean,tempMean - tempErr,'mix_eDown',fig11,'r','--','none')

    
   % r-squared between disc and image
    figure(10); clf;
    subplot(311)
    plot(p_none,p_nat,'go'); hold on;
    plot([0 1],[0 1],'k--')
    xlabel('R2 none'); ylabel('R2 nat');
    subplot(312)
    plot(p_none,p_mix,'ro'); hold on;
    plot([0 1],[0 1],'k--')
    xlabel('R2 none'); ylabel('R2 mix');
    subplot(313)
    plot(p_nat,p_mix,'ko'); hold on;
    plot([0 1],[0 1],'k--')
    xlabel('R2 nat'); ylabel('R2 mix');
    
    [h, p] = ttest(p_none,p_nat);
    [h, p] = ttest(p_none,p_mix);
    [h, p] = ttest(p_nat,p_mix);
    
    
    %Mean NLI comparing surround conditions
    meanNone = mean(meanNLI(:,1)); errNone = std(meanNLI(:,1)) ./ sqrt(length(meanNLI(:,1)));
    meanNat = mean(meanNLI(:,2)); errNat  = std(meanNLI(:,2)) ./ sqrt(length(meanNLI(:,2)));
    meanMix = mean(meanNLI(:,3)); errMix  = std(meanNLI(:,3)) ./ sqrt(length(meanNLI(:,3)));
    
    %No surround vs. natural surround:
    figure; clf;
    fig2=gca;
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','mean NLI none')
    set(get(fig2,'YLabel'),'String','mean NLI nat')
    set(gcf, 'WindowStyle', 'docked')
    addLineToAxis(meanNLI(:,1),meanNLI(:,2),'data',fig2,'g','none','o')
    addLineToAxis(meanNone,meanNat,'mean',fig2,'g','none','s')
    addLineToAxis(meanNone + [-errNone, errNone],[meanNat, meanNat],'errX',fig2,'g','-','none')
    addLineToAxis([meanNone, meanNone],meanNat + [-errNat, errNat],'errY',fig2,'g','-','none')
    addLineToAxis([-0.5 1],[-0.5 1],'unity',fig2,'k','--','none')
    
    %No surround vs. mixed surround:
    figure; clf;
    fig3=gca;
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','mean NLI none')
    set(get(fig3,'YLabel'),'String','mean NLI mixed')
    set(gcf, 'WindowStyle', 'docked')
    addLineToAxis(meanNLI(:,1),meanNLI(:,3),'data',fig3,'r','none','o')
    addLineToAxis(meanNone,meanMix,'mean',fig3,'r','none','s')
    addLineToAxis(meanNone + [-errNone, errNone],[meanMix, meanMix],'errX',fig3,'r','-','none')
    addLineToAxis([meanNone, meanNone],meanMix + [-errMix, errMix],'errY',fig3,'r','-','none')
    addLineToAxis([-0.5 1],[-0.5 1],'unity',fig3,'k','--','none')
    
    %Nat surround vs. mixed surround:
    figure; clf;
    fig4=gca;
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','mean NLI nat')
    set(get(fig4,'YLabel'),'String','mean NLI mixed')
    set(gcf, 'WindowStyle', 'docked')
    addLineToAxis(meanNLI(:,2),meanNLI(:,3),'data',fig4,'k','none','o')
    addLineToAxis(meanNat,meanMix,'mean',fig4,'k','none','s')
    addLineToAxis(meanNat + [-errNat, errNat],[meanMix, meanMix],'errX',fig4,'k','-','none')
    addLineToAxis([meanNat, meanNat],meanMix + [-errMix, errMix],'errY',fig4,'k','-','none')
    addLineToAxis([-0.5 1],[-0.5 1],'unity',fig4,'k','--','none')
    

    % NLI vs mean center intensity
    noBins = 12;
    keepInds = (find(~isnan(allNLIMatrix(:,1))));
    binAndPlotEquallyPopulatedBins(allEqContrast(keepInds),allNLIMatrix(keepInds,1),...
        noBins,fig9,'k','none')    
    
    keepInds = (find(~isnan(allNLIMatrix(:,2))));
    binAndPlotEquallyPopulatedBins(allEqContrast(keepInds),allNLIMatrix(keepInds,2),...
        noBins,fig9,'g','nat')
    
    keepInds = (find(~isnan(allNLIMatrix(:,3))));
    binAndPlotEquallyPopulatedBins(allEqContrast(keepInds),allNLIMatrix(keepInds,3),...
        noBins,fig9,'r','mix')
    

    %paired t-test
    [~, p] = ttest(meanNLI(:,1),meanNLI(:,2))
    [~, p] = ttest(meanNLI(:,1),meanNLI(:,3))
    [~, p] = ttest(meanNLI(:,2),meanNLI(:,3))
    
    
    %Surround conditions cumulative histogram
    edges = linspace(-1,1,40);
    binCtrs = edges(1:end - 1) + mean(diff(edges));
    [nn_none, ~] = histcounts(allNLIMatrix(:,1),edges,'normalization','probability');
    [nn_nat, ~] = histcounts(allNLIMatrix(:,2),edges,'normalization','probability');
    [nn_mix, ~] = histcounts(allNLIMatrix(:,3),edges,'normalization','probability');
    
    %kruskal-wallis test to see if Nat & Shuffled NLI samples come from
    %same distribution
    p = kruskalwallis(allNLIMatrix(:,2:3))
    
    disp('No. patches:')
    disp(size(allNLIMatrix,1));
    disp('No. images:')
    disp(size(meanNLI,1));
    disp('No. cells:')
    disp(pp);
    
    figure; clf;
    fig5=gca;
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','NLI')
    set(get(fig5,'YLabel'),'String','Cumulative prob.')
    set(gcf, 'WindowStyle', 'docked')
    addLineToAxis(binCtrs,cumsum(nn_none),'none',fig5,'k','-','none')
    addLineToAxis(binCtrs,cumsum(nn_nat),'nat',fig5,'g','-','none')
    addLineToAxis(binCtrs,cumsum(nn_mix),'mix',fig5,'r','-','none')

    
    if (exportFigs)
    makeAxisStruct(fig2,'MixSur_meanNLI_1' ,'RFSurroundFigs')
    makeAxisStruct(fig3,'MixSur_meanNLI_2' ,'RFSurroundFigs')
    makeAxisStruct(fig4,'MixSur_meanNLI_3' ,'RFSurroundFigs')
    makeAxisStruct(fig5,'MixSur_NLIcumHist' ,'RFSurroundFigs')
    
    makeAxisStruct(fig6,'MixSur_none' ,'RFSurroundFigs')
    makeAxisStruct(fig7,'MixSur_nat' ,'RFSurroundFigs')
    makeAxisStruct(fig8,'MixSur_mix' ,'RFSurroundFigs')
    
    makeAxisStruct(fig9,'MixSur_vsCen' ,'RFSurroundFigs')
    
    makeAxisStruct(fig11,'MixSur_RvsCen' ,'RFSurroundFigs')
    
    end
    
end