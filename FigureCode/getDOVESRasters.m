function getDOVESRasters(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    
    figColors = pmkmp(8);

    populationNodes = {}; %pop. node is a stimulus
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                '@(list)splitOnRecKeyword(list)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    cell_count = 0;
    figure; clf;
    
    cellData = {};

    cellIDs = {};
    node_count = 0;
    oldCellCount = 0;
    stim_count = 0;
    for pp = 1:length(populationNodes) %for cell in population
        currentNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(currentNode.epochList);
        spikesNode = currentNode.childBySplitValue('extracellular');
        if isempty(spikesNode)
           continue 
        end
        node_count = node_count + 1;
        cellIDs{node_count} = cellInfo.cellID;
        
        cell_count = length(unique(cellIDs));

        if cell_count > oldCellCount
            oldCellCount = cell_count;
            stim_count = 1;
        else
            stim_count = stim_count + 1;
        end

        recType = getRecordingTypeFromEpochList(spikesNode.epochList);
        spikes_centerRes = getMeanResponseTrace(spikesNode.childBySplitValue('Center').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        spikes_surroundRes = getMeanResponseTrace(spikesNode.childBySplitValue('Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        spikes_centerSurroundRes = getMeanResponseTrace(spikesNode.childBySplitValue('Center-Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
        
        % save out spike binary data...
        new_cell_struct = struct;
        new_cell_struct.cellID = cellInfo.cellID;
        new_cell_struct.cellType = cellInfo.cellType;
        new_cell_struct.stimulusID = currentNode.splitValue;
        new_cell_struct.Center = spikes_centerRes.binary;
        new_cell_struct.Surround = spikes_surroundRes.binary;
        new_cell_struct.CenterSurround = spikes_centerSurroundRes.binary;
        cellData{node_count} = new_cell_struct;
        
        %timing stuff:
        timingEpoch = currentNode.epochList.firstValue;
        frameRate = timingEpoch.protocolSettings('background:Microdisplay Stage@localhost:monitorRefreshRate');
        if isempty(frameRate)
            frameRate = timingEpoch.protocolSettings('background:Microdisplay_Stage@localhost:monitorRefreshRate');
        end
        FMdata = (riekesuite.getResponseVector(timingEpoch,'Frame Monitor'))';
        [measuredFrameTimes, ~] = getFrameTiming(FMdata,0);

        sampleRate = currentNode.epochList.firstValue.protocolSettings('sampleRate');
        preTime = currentNode.epochList.firstValue.protocolSettings('preTime') / 1e3; %sec
        stimTime = currentNode.epochList.firstValue.protocolSettings('stimTime') / 1e3; %sec
        preFrames = preTime .* frameRate;
        firstStimFrameFlip = measuredFrameTimes(preFrames + 1);
        
        %Stimulus stuff. To get fixation periods...
        stimSet = currentNode.epochList.firstValue.protocolSettings('currentStimSet');
        imageSet = 'NaturalImages';
        resourcesDir = '~/Dropbox/RiekeLab/Analysis/MATLAB/turner-package/resources/';
        load([resourcesDir, stimSet]);
        eyeX = FEMdata(currentNode.splitValue).eyeX;
        eyeY = FEMdata(currentNode.splitValue).eyeY;
        dx = abs(diff(eyeX));
        dy = abs(diff(eyeY));
        instVel = sqrt(dx.^2 + dy.^2); %instantaneous velocity
        
        [peaks,peakInd] = getPeaks(instVel,1);
        threshInd = peaks > 10; %only changes greater than 10 amin / unit time
        peakInd = peakInd(threshInd);
        
        %10 sample points at 200 Hz -> 50 msec apart. Filters out quick
        %little "compound" saccadic movements
        refractoryPeriod = 10;
        filtInds = peakInd(1);
        ct = 1;
        for ii = 2:length(peakInd)
            if peakInd(ii) - filtInds(ct) < refractoryPeriod
            else
                ct = ct + 1;
                filtInds = cat(1,filtInds,peakInd(ii));
            end
        end
        filtInds = [0 filtInds']; %start at 1 and shift back due to diff above
        filtInds(end + 1) = length(eyeX); %end of last fixation
        saccadeTimes = round((filtInds ./ 200) .* sampleRate + firstStimFrameFlip); %datapoints

%         figure; clf; fig1=gca; initFig(fig1,'Time (s)','Trial')
%         subplot(10,3,(cell_count-1)*3 + stim_count)
%         fh = gca;
%         addRastersToFigure(spikes_centerSurroundRes.binary,fh)

%         addRastersToFigure(spikes_centerRes.binary,fig1)                    
%         addRastersToFigure(spikes_surroundRes.binary,fig12)                    

        %natural image with overlayed eye movements
%         imageName = FEMdata(currentNode.splitValue).ImageName;
%         fileId=fopen([resourcesDir, imageSet, '/', imageName],'rb','ieee-be');
%         img = fread(fileId, [1536,1024], 'uint16');
%         img = double(img);
%         img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
%         fh = figure(21); clf;
%         imagesc(img'); colormap(gray); axis image; axis off; hold on;
%         plot(FEMdata(currentNode.splitValue).eyeX,FEMdata(currentNode.splitValue).eyeY,'r')
%         brighten(0.6) %brighten for display purposes
%         drawnow;
%         figID = 'egNatImage_DOVES';
%         print(fh,[figDir,figID],'-depsc')
   
    end %end for population
    
    save('DovesRastersData.mat','cellData')
end