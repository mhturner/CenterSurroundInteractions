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
    
    for pp = 1:length(populationNodes) % for cell in pop
        cellNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(cellNode.epochList);
        recType = getRecordingTypeFromEpochList(cellNode.epochList);
        
        ImageResponseMatrix = []; % rows = center image, columns = surround image ([none, natural, mixed1, mixed2,...])
        DiscResponseMatrix = [];
        NLIresultsMatrix = []; 
        for imageInd = 1:cellNode.children.length % for image ID
            imageNode = cellNode.children(imageInd);
            for patchInd = 1:imageNode.children.length % for center patch
                patchNode = imageNode.children(patchInd);
                centerPatchLocation = convertJavaArrayList(patchNode.epochList.firstValue.protocolSettings('centerPatchLocation'));
                mixedSurroundCount = 0;
                for surroundInd = 1:patchNode.children.length % for surround patch
                    surroundNode = patchNode.children(surroundInd);
                    currentSurroundLocation = convertJavaArrayList(surroundNode.epochList.firstValue.protocolSettings('currentSurroundLocation'));
                    %where to put results in NLIresultsMatrix:
                    if isequal(currentSurroundLocation, [0, 0]) %no surround
                        putInd = 1;
                    elseif isequal(currentSurroundLocation, centerPatchLocation) %natural surround
                        putInd = 2;
                    else %mixed surround
                        mixedSurroundCount = mixedSurroundCount + 1;
                        putInd = 2 + mixedSurroundCount;
                    end
                    
                    %get image and disc responses:
                    imageResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('image').epochList,recType);
                    discResp = getResponseAmplitudeStats(surroundNode.childBySplitValue('intensity').epochList,recType);
                    %compute NLI:
                    newDiff = (imageResp.(metric).mean - discResp.(metric).mean);
                    newNLIvalue =  newDiff ./ ...
                           (abs(imageResp.(metric).mean) + abs(discResp.(metric).mean));
                       
                    ImageResponseMatrix(patchInd,putInd) = imageResp.(metric).mean;
                    DiscResponseMatrix(patchInd,putInd) = discResp.(metric).mean;
                    NLIresultsMatrix(patchInd,putInd) = newNLIvalue;
                    
                end % for surround patch
                
                
            end % for center patch
            
        end % for image ID
    end % for cell in pop
    
    NLI_rand = mean(NLIresultsMatrix(:,3:end),2);
    NLI_nat = NLIresultsMatrix(:,2);
    NLI_no = NLIresultsMatrix(:,1);
    
    figure(9); clf;
    subplot(221); plot(NLI_no,NLI_nat,'go'); hold on; plot([0 1],[0 1],'k--')
    subplot(222); plot(NLI_no,NLI_rand,'ro'); hold on; plot([0 1],[0 1],'k--')
    subplot(223); plot(NLI_nat,NLI_rand,'ko'); hold on; plot([0 1],[0 1],'k--')

    figure(10); clf; 
    hold on;

    noPatches = size(NLIresultsMatrix,1);
    noRandSurrounds = size(NLIresultsMatrix,2) - 2;
    colors = pmkmp(noPatches);
    for patchInd = 1:noPatches
        tempNatSurround = NLIresultsMatrix(patchInd,1);
        plot(repmat(tempNatSurround,1,noRandSurrounds),NLIresultsMatrix(patchInd,3:end),...
            'Color',colors(patchInd,:),'Marker','o','LineStyle','none')
        
        plot(NLIresultsMatrix(patchInd,1),NLIresultsMatrix(patchInd,2),...
            'Color',colors(patchInd,:),'Marker','x','LineStyle','none');

    end

    hold on; plot([0 1],[0 1],'k--')
    xlabel('NLI, no surround'); ylabel('NLI, with surround');
    
    noBins = 5;
    figure(11); clf; hold on;
    [nn, ctr] = histcounts(NLIresultsMatrix(:,1),noBins);
    plot(ctr(2:end),nn./sum(nn),'k');
    
    [nn, ctr] = histcounts(NLIresultsMatrix(:,2),noBins);
    plot(ctr(2:end),nn./sum(nn),'g');
    
    temp = NLIresultsMatrix(:,3:end);
    temp = temp(:);
    [nn, ctr] = histcounts(temp,noBins);
    plot(ctr(2:end),nn./sum(nn),'r');
    
    
end