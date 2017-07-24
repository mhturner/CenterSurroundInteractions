function doLEDMixedSurroundAnalysis(node,varargin)
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
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(imageName)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    meanNLI = [];
    allNLIMatrix = [];
    
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
        for imageInd = 1:cellNode.children.length % for image ID
            imageCt = imageCt + 1;
            imageNode = cellNode.children(imageInd);
            for patchInd = 1:imageNode.children.length % for center patch
                patchNode = imageNode.children(patchInd);                
                
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
            
            %compute NLI matrix for this image
            % rows = center patch, 3 columns = surround condition ([none, natural, mixed])
            diffMatrix = ImageResponseMatrix - DiscResponseMatrix;
            NLIresultsMatrix = diffMatrix ./ (abs(ImageResponseMatrix) + abs(DiscResponseMatrix));
            
            allNLIMatrix = cat(1,allNLIMatrix,NLIresultsMatrix);
            
            meanNLI(imageCt,:) = nanmean(NLIresultsMatrix,1);

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
    
   
    figure(4); clf;
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
    
    
    figure(5); clf;
    subplot(311)
    plot(meanNLI(:,1),meanNLI(:,2),'go'); hold on;
    plot([0 1],[0 1],'k--')
    xlabel('mean NLI none'); ylabel('mean NLI nat');
    subplot(312)
    plot(meanNLI(:,1),meanNLI(:,3),'ro'); hold on;
    plot([0 1],[0 1],'k--')
    xlabel('mean NLI none'); ylabel('mean NLI mix');
    subplot(313)
    plot(meanNLI(:,2),meanNLI(:,3),'ko'); hold on;
    plot([0 1],[0 1],'k--')
    xlabel('mean NLI nat'); ylabel('mean NLI mix');
    
    [h, p] = ttest(meanNLI(:,1),meanNLI(:,2))
    [h, p] = ttest(meanNLI(:,1),meanNLI(:,3))
    [h, p] = ttest(meanNLI(:,2),meanNLI(:,3))
    
    
    
    edges = linspace(-1,1,40);
    binCtrs = edges(1:end - 1) + mean(diff(edges));
    [nn_none, ~] = histcounts(allNLIMatrix(:,1),edges,'normalization','probability');
    [nn_nat, ~] = histcounts(allNLIMatrix(:,2),edges,'normalization','probability');
    [nn_mix, ~] = histcounts(allNLIMatrix(:,3),edges,'normalization','probability');
    
    figure(6); clf; hold on;
    plot(binCtrs,cumsum(nn_none),'k')
    plot(binCtrs,cumsum(nn_nat),'g')
    plot(binCtrs,cumsum(nn_mix),'r')

    
end