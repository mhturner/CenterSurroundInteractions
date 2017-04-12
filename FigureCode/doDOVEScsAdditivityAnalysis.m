function doDOVEScsAdditivityAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'exportFigs',true,@islogical);
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    exportFigs = ip.Results.exportFigs;
    
    figColors = pmkmp(8);
    figDir = '~/Documents/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs
    
    figure; clf;
    fig1=gca; % Response traces, c and s
    set(fig1,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig1,'XLabel'),'String','Time (s)')
    set(get(fig1,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf;
    fig2=gca; % Response traces, linear sum and cs
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Time (s)')
    set(get(fig2,'YLabel'),'String','Response')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf;
    fig3=gca; % Eye movement positions over time
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Time (s)')
    set(get(fig3,'YLabel'),'String','Location (arcmin)')
    set(gcf, 'WindowStyle', 'docked')

    figure; clf;
    fig11=gca; % RASTER - center alone (spike rec only)
    set(fig11,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig11,'XLabel'),'String','Time (s)')
    set(get(fig11,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf;
    fig12=gca; % RASTER - surround alone (spike rec only)
    set(fig12,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig12,'XLabel'),'String','Time (s)')
    set(get(fig12,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf;
    fig13=gca; % RASTER - center + surround (spike rec only)
    set(fig13,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig13,'XLabel'),'String','Time (s)')
    set(get(fig13,'YLabel'),'String','Trial')
    set(gcf, 'WindowStyle', 'docked')
    
    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(stimulusIndex)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end

    for pp = 1:length(populationNodes) %for cell in population
        currentNode = populationNodes{pp};
        recType = getRecordingTypeFromEpochList(currentNode.epochList);
        cellInfo = getCellInfoFromEpochList(currentNode.epochList);
        for ss = 1:currentNode.children.length %for stim
            stimNode = currentNode.children(ss);
            centerRes = getMeanResponseTrace(stimNode.childBySplitValue('Center').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
            surroundRes = getMeanResponseTrace(stimNode.childBySplitValue('Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
            centerSurroundRes = getMeanResponseTrace(stimNode.childBySplitValue('Center-Surround').epochList,recType,'PSTHsigma',10,'attachSpikeBinary',true);
            
            if stimNode.custom.get('isExample')
                addLineToAxis(centerRes.timeVector,centerRes.mean,...
                    'center',fig1,figColors(1,:),'-','none')
                addLineToAxis(surroundRes.timeVector,surroundRes.mean,...
                    'surround',fig1,figColors(4,:),'-','none')
                addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(stimNode.splitValue)],fig1,'k','none','none')
                
                addLineToAxis(centerRes.timeVector,centerRes.mean + surroundRes.mean,...
                    'linSum',fig2,[0.7 0.7 0.7],'-','none')
                addLineToAxis(centerSurroundRes.timeVector,centerSurroundRes.mean,...
                    'measuredCS',fig2,'k','-','none')
                addLineToAxis(0,0,[cellInfo.cellID,'_', num2str(stimNode.splitValue)],fig2,'k','none','none')
                if strcmp(recType,'extracellular') %do raster figs
                    addRastersToFigure(centerRes.binary,fig11)                    
                    addRastersToFigure(surroundRes.binary,fig12)                    
                    addRastersToFigure(centerSurroundRes.binary,fig13)
                end
                
                stimSet = stimNode.epochList.firstValue.protocolSettings('currentStimSet');
                imageSet = 'NaturalImages';
                resourcesDir = '~/Documents/MATLAB/turner-package/resources/';
                load([resourcesDir, stimSet]);
                
                %natural image with overlayed eye movements
                imageName = FEMdata(stimNode.splitValue).ImageName;
                fileId=fopen([resourcesDir, imageSet, '/', imageName],'rb','ieee-be');
                img = fread(fileId, [1536,1024], 'uint16');
                img = double(img);
                img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level
                fh = figure(21); clf;
                imagesc(img'); colormap(gray); axis image; axis off; hold on;
                plot(FEMdata(stimNode.splitValue).eyeX,FEMdata(stimNode.splitValue).eyeY,'r')
                drawnow;
                figID = 'egNatImage_DOVES';
                print(fh,[figDir,figID],'-depsc')
                
                %eye position over time
                preTime = stimNode.epochList.firstValue.protocolSettings('preTime') / 1e3; %sec
                stimTime = stimNode.epochList.firstValue.protocolSettings('stimTime') / 1e3; %sec
                tailTime = stimNode.epochList.firstValue.protocolSettings('tailTime') / 1e3; %sec
                totalTime = preTime + stimTime + tailTime; %sec
                
                totalOnPoints = stimTime * 200;
                movementPoints = length(FEMdata(stimNode.splitValue).eyeX); %at 200 Hz
                xTrace = FEMdata(stimNode.splitValue).eyeX;
                yTrace = FEMdata(stimNode.splitValue).eyeY;
                if movementPoints < totalOnPoints
                   xTrace((movementPoints + 1):totalOnPoints) = FEMdata(stimNode.splitValue).eyeX(end);
                   yTrace((movementPoints + 1):totalOnPoints) = FEMdata(stimNode.splitValue).eyeY(end);
                end
                
                onTimeVec = preTime + (1:totalOnPoints) ./ 200;
                addLineToAxis(onTimeVec,xTrace,...
                    'xTrace',fig3,figColors(2,:),'-','none')
                addLineToAxis(onTimeVec,yTrace,...
                    'yTrace',fig3,figColors(6,:),'-','none')
                addLineToAxis([onTimeVec(1) onTimeVec(1)],[0 1200],...
                    'imageON',fig3,[0.7 0.7 0.7],'--','none')
                addLineToAxis([onTimeVec(end) onTimeVec(end)],[0 1200],...
                    'imageOFF',fig3,[0.7 0.7 0.7],'--','none')
            end %end if example
            
        end %end for stim
        
    end %end for population
    
    if (exportFigs)
        figID = ['DOVEScs_ind_',recType];
        makeAxisStruct(fig1,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_CS_',recType];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')
        
        figID = ['DOVEScs_em'];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')
        
        if strcmp(recType,'extracellular') %raster plots
            figID = 'DOVEScs_rasterC';
            makeAxisStruct(fig11,figID ,'RFSurroundFigs')
            
            figID = 'DOVEScs_rasterS';
            makeAxisStruct(fig12,figID ,'RFSurroundFigs')
            
            figID = 'DOVEScs_rasterCS';
            makeAxisStruct(fig13,figID ,'RFSurroundFigs')
        end
    end
    
end