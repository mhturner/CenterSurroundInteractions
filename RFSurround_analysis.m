%% RESET. NEW TREE.
clear list; clear all; clear java; close all; clc; %#ok<CLJAVA,CLALL>
loader = edu.washington.rieke.Analysis.getEntityLoader();
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();
dataFolder = '/Users/maxturner/CurrentData/RFSurround/'; 
saveFileDirectory = '~/Documents/MATLAB/Analysis/Projects/RFSurround/SavedTreeFlags/';
import auimodel.*
import vuidocument.*
cd('~/Documents/MATLAB/Analysis/Projects/RFSurround/')
%% EXPANDING SPOTS: tree
% saveFileID = 'ES_HorizontalCellNegative';

list = loader.loadEpochList([dataFolder,'ExpandingSpots.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,centerSplit_java,...
    'protocolSettings(spotIntensity)','protocolSettings(currentSpotSize)'});

gui = epochTreeGUI(tree);
if exist('saveFileID','var')
    setFlaggedAndExampleNodes(gui, tree, saveFileDirectory, saveFileID)
end
%% EXPANDING SPOTS: flag population and set example(s)
%   Select cell type as parentNode
%       Flag spotContrast node of cells in population    

% save flagged and example-d tree nodes for the future
saveFileID = 'ES_HorizontalCellNegative';
getFlaggedAndExampleNodes(tree, saveFileDirectory, saveFileID);
%% EXPANDING SPOTS: do analysis & make figs
%       Select cell type as parentNode
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doAreaSummationAnalysis(parentNode,'metric','integrated',...
    'amplitudeMultiplier',1);

%% CONTRAST-REVERSING GRATINGS: tree
saveFileID = 'CRGs_exc';

list = loader.loadEpochList([dataFolder,'CRG-noDrug.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

maskSplit = @(list)splitOnRadiusOrDiameter(list,'mask');
maskSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, maskSplit);

apertureSplit = @(list)splitOnRadiusOrDiameter(list,'aperture');
apertureSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, apertureSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,centerSplit_java,...
    maskSplit_java,...
    apertureSplit_java,...
    'protocolSettings(currentBarWidth)'});

gui = epochTreeGUI(tree);
if exist('saveFileID','var')
    setFlaggedAndExampleNodes(gui, tree, saveFileDirectory, saveFileID)
end
%% CONTRAST-REVERSING GRATINGS: get population
%   Select cell type as parentNode
%       Flag apertureDiamter node of cells in population    

% save flagged and example-d tree nodes for the future
saveFileID = 'CRGs_exc';
getFlaggedAndExampleNodes(tree, saveFileDirectory, saveFileID);
%% CONTRAST-REVERSING GRATINGS: do analysis & make figs
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doContrastReversingGratingsAnalysis(parentNode,...
    'normalizePopulationF2',true,'noBins',10);

%% LINEAR EQUIVALENT ANNULUS: tree
list = loader.loadEpochList([dataFolder,'LinearEquivalentAnnulus.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

patchSamplingSplit = @(list)splitOnPatchSampling_NatImage(list);
patchSamplingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, patchSamplingSplit);

patchContrastSplit = @(list)splitOnPatchContrast_NatImage(list);
patchContrastSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, patchContrastSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'currentPatchLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',recordingSplit_java,...
    patchSamplingSplit_java,...
    patchContrastSplit_java,...
    'protocolSettings(rfSigmaSurround)',...
    'protocolSettings(annulusInnerDiameter)',...
    'protocolSettings(annulusOuterDiameter)',...
    'protocolSettings(linearIntegrationFunction)',...
    'protocolSettings(centerSpotContrast)',...
    'protocolSettings(imageName)',...
    locationSplit_java,...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% LINEAR EQUIVALENT ANNULUS: select example nodes and make figs
%   Select cell type as parentNode
%       Flag linearIntegrationFunction node of cells in population    
clc;
% parentNode = gui.getSelectedEpochTreeNodes{1};
doLinearEquivalentAnalysis(parentNode,...
    'metric','integrated','stimType','A',...
    'subunitSigma',25,...
    'centerSigma',150);

%% LINEAR EQUIVALENT DISC: tree
list = loader.loadEpochList([dataFolder,'LinearEquivalentDisc.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

patchSamplingSplit = @(list)splitOnPatchSampling_NatImage(list);
patchSamplingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, patchSamplingSplit);

patchContrastSplit = @(list)splitOnPatchContrast_NatImage(list);
patchContrastSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, patchContrastSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'currentPatchLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',recordingSplit_java,...
    patchSamplingSplit_java,...
    patchContrastSplit_java,...
    'protocolSettings(apertureDiameter)',...
    'protocolSettings(linearIntegrationFunction)',...
    'protocolSettings(imageName)',...
    locationSplit_java,...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);
%% LINEAR EQUIVALENT DISC: get population
%   Select cell type as parentNode
%       Flag linearIntegrationFunction node of cells in population    
parentNode = gui.getSelectedEpochTreeNodes{1};
%% LINEAR EQUIVALENT DISC: select example nodes and make figs
clc;
doLinearEquivalentAnalysis(parentNode,...
    'metric','integrated','stimType','D');

%% NatImageCSAdditivity tree

list = loader.loadEpochList([dataFolder,'NatImageCSAdditivity.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'currentPatchLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,centerSplit_java,...
    'protocolSettings(patchContrast)',...
    'protocolSettings(patchSampling)',...
    'protocolSettings(imageName)',...
    locationSplit_java,...
    'protocolSettings(currentStimulus)'});

gui = epochTreeGUI(tree);
%% NatImageCSAdditivity simple plots
%select image name node
parentNode = gui.getSelectedEpochTreeNodes{1};

responses = zeros(1,3); %c, s, cs
ct = 0;
for pp = 1:parentNode.children.length
    imageNode = parentNode.children(pp);
    for ll = 1:parentNode.children(pp).children.length
        locationNode = parentNode.children(pp).children(ll);
        if locationNode.children.length < 3
            continue
        end
        ct = ct + 1;
        tempCenter = getResponseAmplitudeStats(locationNode.childBySplitValue('Center').epochList,'extracellular');
        responses(ct,1) = tempCenter.integrated.mean;

        tempSurround = getResponseAmplitudeStats(locationNode.childBySplitValue('Surround').epochList,'extracellular');
        responses(ct,2) = tempSurround.integrated.mean;

        tempCS = getResponseAmplitudeStats(locationNode.childBySplitValue('Center-Surround').epochList,'extracellular');
        responses(ct,3) = tempCS.integrated.mean;

    end
end

%%
figure(2); clf; hold on;
[a, s_c] = getActivityRatio(responses(:,1));
[n, edges] = histcounts(responses(:,1),5);
centers = diff(edges)./2 + edges(1:end-1);
plot(centers, n, 'b')

[a, s_cs] = getActivityRatio(responses(:,3));
[n, edges] = histcounts(responses(:,3),5);
centers = diff(edges)./2 + edges(1:end-1);
plot(centers, n, 'k')


    
%%
figure(2); clf;
subplot(221); hold on;
plot(responses(:,1) + responses(:,2),responses(:,3),'ko')
plot([min(responses(:)) max(responses(:))],[min(responses(:)) max(responses(:))],'k--')
xlabel('R(C) + R(S)'); ylabel('R(C+S)')

subplot(222); hold on;
plot(responses(:,1),responses(:,2),'ko')
plot([min(responses(:)) max(responses(:))],[min(responses(:)) max(responses(:))],'k--')
xlabel('R(C)'); ylabel('R(S)')

subplot(223); hold on;
plot(responses(:,1),responses(:,3),'ko')
plot([min(responses(:)) max(responses(:))],[min(responses(:)) max(responses(:))],'k--')
xlabel('R(C)'); ylabel('R(C+S)')


%% F1F2 CONTRAST: tree
originalList = loader.loadEpochList([dataFolder,'contrastF1F1&SFC.mat'],dataFolder);
saveFileID = 'F2F1_exc';

for ii = 1:2
    if ii == 1
        list = originalList;
    elseif ii == 2
        list = newEL;
    end
    recordingSplit = @(list)splitOnRecKeyword(list);
    recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

    cellTypeSplit = @(list)splitOnCellType(list);
    cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

    centerSizeSplit = @(list)splitOnF1F2CenterSize(list);
    centerSizeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSizeSplit);

    contrastSplit = @(list)splitOnF1F2Contrast(list);
    contrastSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, contrastSplit);

    phaseSplit = @(list)splitOnF1F2Phase(list);
    phaseSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, phaseSplit);

    tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
        recordingSplit_java,...
        'protocolSettings(maskDiameter)',...
        centerSizeSplit_java,...
        phaseSplit_java,...
        contrastSplit_java});

    if ii == 1
        newEL = edu.washington.rieke.Analysis.getEpochListFactory().create;
        for nn = 1:tree.descendentsDepthFirst.length
            if tree.descendentsDepthFirst(nn).splitKey == contrastSplit_java
                if tree.descendentsDepthFirst(nn).children.length > 1
                    for ee = 1:tree.descendentsDepthFirst(nn).epochList.length
                        newEL.append(tree.descendentsDepthFirst(nn).epochList.elements(ee));
                    end
                end
            end
        end
    end
    
end
gui = epochTreeGUI(tree);
if exist('saveFileID','var')
    setFlaggedAndExampleNodes(gui, tree, saveFileDirectory, saveFileID)
end
%% F1F2 CONTRAST: get population
%   Select cell type as parentNode
%       Flag center size nodes. For center-surround cells, flag both center
%       sizes

% save flagged and example-d tree nodes for the future
saveFileID = 'F2F1_exc';
getFlaggedAndExampleNodes(tree, saveFileDirectory, saveFileID);

%% F1F2 CONTRAST: select example nodes and make figs
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doContrastF1F2Analysis(parentNode,phaseSplit_java,'noBins',6);

%% HEPES EXPANDING SPOTS: tree

list_full = loader.loadEpochList([dataFolder,'ExpandingSpots.mat'],dataFolder);

targetGroups = {'NaCl BA control','HEPES'};
list = filterEpochListByEpochGroups(list_full,targetGroups,'Include');

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,centerSplit_java,'protocolSettings(spotIntensity)',...
    'protocolSettings(epochGroup:label)',...
    'protocolSettings(currentSpotSize)'});
gui = epochTreeGUI(tree);


%% HEPES EXPANDING SPOTS: do analysis & make figs
%       Select cell type as parentNode
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doAreaSummationAnalysis(parentNode,'metric','integrated',...
    'figureID','OFFparHEPES1');

%% HEPES CONTRAST-REVERSING GRATINGS: tree

list_full = loader.loadEpochList([dataFolder,'CRG.mat'],dataFolder);

targetGroups = {'NaCl BA control','HEPES'};
list = filterEpochListByEpochGroups(list_full,targetGroups,'Include');

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

maskSplit = @(list)splitOnRadiusOrDiameter(list,'mask');
maskSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, maskSplit);

apertureSplit = @(list)splitOnRadiusOrDiameter(list,'aperture');
apertureSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, apertureSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,centerSplit_java,...
    maskSplit_java,...
    apertureSplit_java,'protocolSettings(epochGroup:label)',...
    'protocolSettings(currentBarWidth)'});

gui = epochTreeGUI(tree);

%% CONTRAST-REVERSING GRATINGS: do analysis & make figs
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doContrastReversingGratingsAnalysis(parentNode,...
    'normalizePopulationF2',true,'noBins',10,...
    'figureID','ONparControl');


%% TTX CRGs: tree
list = loader.loadEpochList([dataFolder,'CRG.mat'],dataFolder);
recordingSplit = @(list)splitOnKeywords(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

tree = riekesuite.analysis.buildTree(list, {'protocolSettings(source:type)','cell.label',recordingSplit_java,...
    'protocolSettings(maskDiameter)',...
    'protocolSettings(apertureDiameter)',...
    'protocolSettings(currentBarWidth)'});

%check that some parameters are consistent within cells...
constantSettings = {'backgroundIntensity','preTime','stimTime','tailTime','contrast','temporalFrequency'};
for cellTypeIndex = 1:tree.children.length
    for cellNameIndex = 1:tree.children(cellTypeIndex).children.length
        newEL = tree.children(cellTypeIndex).children(cellNameIndex).epochList;
        [outEpochList , ~] = makeUniformEpochList(newEL,constantSettings,[]);
    end
end

gui = epochTreeGUI(tree);

%% TTX CRGs: example cell plotting
%   Select recording type (exc, inh, etc) as parentNode
clc
parentNode = gui.getSelectedEpochTreeNodes{1};
cellInfo = getCellInfoFromEpochList(parentNode.epochList);

%get bar widths for color code
bw = zeros(1,parentNode.epochList.length);
for ee = 1:parentNode.epochList.length
   bw(ee) = parentNode.epochList.elements(ee).protocolSettings('currentBarWidth');
end
uniqueBarWidths = unique(bw);

drug = []; wash = []; preDrug = [];
for cc = 1:parentNode.children.length
    if parentNode.children(cc).custom.get('isSelected')
        for ii = 1:parentNode.children(cc).children.length
            if parentNode.children(cc).children(ii).splitValue > 0
                node = parentNode.children(cc).children(ii).children(1); % SURROUND: MASK > 0
                break
            end
        end
    else
        continue
    end

    tempWid = nan(1,node.children.length); tempF2 = nan(1,node.children.length);
    for CurNode = 1:node.children.length
        currentBarWidth = node.children(CurNode).splitValue;
        recType = getRecordingTypeFromEpochList(node.children(CurNode).epochList);
        responseTrace = getCycleAverageResponse(node.children(CurNode).epochList,recType,preTime,temporalFrequency);

        stats = getF1F2statistics(node.children(CurNode).epochList,recType,preTime,temporalFrequency);
        tempF2(CurNode) = stats.meanF2;
        tempWid(CurNode) = currentBarWidth;
    end

    conditionKeywords = node.parent.parent.splitValue;
    if ~isempty(strfind(conditionKeywords,'nM'))
        drug.barWidth = tempWid; drug.F2 = tempF2;
    elseif ~isempty(strfind(conditionKeywords,'wash'))
        wash.barWidth = tempWid; wash.F2 = tempF2;
    else
        preDrug.barWidth = tempWid; preDrug.F2 = tempF2;
    end
end

figure; clf;
fig1=gca;
set(fig1,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Bar width (um)')
set(get(fig1,'YLabel'),'String','F2 amplitude (norm)')

addLineToAxis(preDrug.barWidth,preDrug.F2 ./ max(preDrug.F2),'preDrug',fig1,'k','-','o')
addLineToAxis(drug.barWidth,drug.F2 ./ max(drug.F2),'drug',fig1,'r','-','o')
addLineToAxis(wash.barWidth,wash.F2 ./ max(wash.F2),'wash',fig1,'b','-','o')

makeAxisStruct(fig1,['TtxF2Sur_',cellInfo.cellType,'_',recType] ,'RFSurroundFigs')

%% HORIZONTAL CELL F2 AT LIGHT LEVEL & VREST: tree
list = loader.loadEpochList([dataFolder,'SplitFieldCentering.mat'],dataFolder);
cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

tempTree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,recordingSplit_java});
list = tempTree.childBySplitValue('horizontal').childBySplitValue('iClamp').epochList;
constrainedSettings.maskDiameter = 0;
constrainedSettings.contrast = 0.9;
[list, outParams] = makeUniformEpochList(list,[],constrainedSettings);


centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

monLevelSplit = @(list)splitOnOLEDLevel(list);
monLevelSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, monLevelSplit);

holdingCurrentSplit = @(list)splitOnHoldingSignal(list);
holdingCurrentSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, holdingCurrentSplit);


tree = riekesuite.analysis.buildTree(list, {'cell.label',...
    centerSplit_java,...
    'protocolSettings(spotDiameter)',...
    monLevelSplit_java,...
    holdingCurrentSplit_java,...
    'protocolSettings(splitField)'});

%check that some parameters are consistent within cells...
constantSettings = {'backgroundIntensity','preTime','stimTime','tailTime'};
for cellTypeIndex = 1:tree.children.length
    for cellNameIndex = 1:tree.children(cellTypeIndex).children.length
        newEL = tree.children(cellTypeIndex).children(cellNameIndex).epochList;
        [outEpochList, outParams] = makeUniformEpochList(newEL,constantSettings,[]);
    end
end

gui = epochTreeGUI(tree);

%% HORIZONTAL CELL F2 AT LIGHT LEVEL & VREST: example cell plotting
%select spot diameter node
clc
parentNode = gui.getSelectedEpochTreeNodes{1};
recType = getRecordingTypeFromEpochList(parentNode.epochList);

% % % First: no holding current, just changing background
medNode = parentNode.childBySplitValue('medium').childBySplitValue(0);
highNode = parentNode.childBySplitValue('high').childBySplitValue(0);

figure; clf;
fig1=gca;
set(fig1,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Time (sec)')
set(get(fig1,'YLabel'),'String','Vm (mV)')

stats_F1 = getF1F2statistics(medNode.children(1).epochList,recType);
stats_F2 = getF1F2statistics(medNode.children(2).epochList,recType);
F2F1Ratio_med = stats_F2.meanF2 / stats_F1.meanF1;
temp_F1 = getMeanResponseTrace(medNode.children(1).epochList,recType);
temp_F2 = getMeanResponseTrace(medNode.children(2).epochList,recType);
addLineToAxis(temp_F1.timeVector,temp_F1.mean,'mediumF1',fig1,'k','-','none')
addLineToAxis(temp_F2.timeVector,temp_F2.mean,'mediumF2',fig1,'r','-','none')


stats_F1 = getF1F2statistics(highNode.children(1).epochList,recType);
stats_F2 = getF1F2statistics(highNode.children(2).epochList,recType);
F2F1Ratio_high = stats_F2.meanF2 / stats_F1.meanF1;
temp_F1 = getMeanResponseTrace(highNode.children(1).epochList,recType);
temp_F2 = getMeanResponseTrace(highNode.children(2).epochList,recType);
addLineToAxis(temp_F1.timeVector,temp_F1.mean,'highF1',fig1,'k','-','none')
addLineToAxis(temp_F2.timeVector,temp_F2.mean,'highF2',fig1,'r','-','none')

str = ['F2/F1, 9,000 = ',num2str(F2F1Ratio_high),'; 500 = ',num2str(F2F1Ratio_med)];
title(str)

% makeAxisStruct(fig1,'MonLevelF1F2' ,'RFSurroundFigs')

%% % % Second: mon high, injecting current:
parentNode = gui.getSelectedEpochTreeNodes{1};
recType = getRecordingTypeFromEpochList(parentNode.epochList);
highNode = parentNode.childBySplitValue('high');


figure; clf;
fig2=gca;
set(fig2,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig2,'XLabel'),'String','Time (sec)')
set(get(fig2,'YLabel'),'String','Vm (mV)')
colors = pmkmp(highNode.children.length);
F1s = nan(1,highNode.children.length);
F2s = nan(1,highNode.children.length);
baselines = nan(1,highNode.children.length);
currents = nan(1,highNode.children.length);
for cc = 1:highNode.children.length
    stats_F1 = getF1F2statistics(highNode.children(cc).children(1).epochList,recType);
    stats_F2 = getF1F2statistics(highNode.children(cc).children(2).epochList,recType);
    
    F1s(cc) = stats_F1.meanF1;
    F2s(cc) = stats_F2.meanF2;

    temp_F1 = getMeanResponseTrace(highNode.children(cc).children(1).epochList,recType);
    temp_F2 = getMeanResponseTrace(highNode.children(cc).children(2).epochList,recType);
    
    baselines(cc) = temp_F1.baseline;
    currents(cc) = highNode.children(cc).splitValue;
    
    addLineToAxis(temp_F1.timeVector,temp_F1.mean,'highF1',fig2,colors(cc,:),'-','none')
    addLineToAxis(temp_F2.timeVector,temp_F2.mean,'highF2',fig2,colors(cc,:),'-','none')
end

figure; clf;
fig3=gca;
set(fig3,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig3,'XLabel'),'String','Vrest (mV)')
set(get(fig3,'YLabel'),'String','Amplitude (mV)')

addLineToAxis(baselines,F1s,'highF1',fig3,'k','-','o')
addLineToAxis(baselines,F2s,'highF1',fig3,'r','-','o')
legend('F1','F2')

figure; clf;
fig4=gca;
set(fig4,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig4,'XLabel'),'String','Vrest (mV)')
set(get(fig4,'YLabel'),'String','F2:F1')
addLineToAxis(baselines,F2s ./ F1s,'highF1',fig4,'k','-','o')

%% HORIZONTAL CELL CURRENT INJECTION: tree
list = loader.loadEpochList([dataFolder,'SinusoidalCurrentInjection.mat'],dataFolder);


monLevelSplit = @(list)splitOnOLEDLevel(list);
monLevelSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, monLevelSplit);

tree = riekesuite.analysis.buildTree(list, {'cell.label',...
    monLevelSplit_java,...
    'protocolSettings(stimuli:Amp1:mean)',...
    'protocolSettings(stimuli:Amp1:amplitude)'});
gui = epochTreeGUI(tree);



%%
%select cell node
clc
parentNode = gui.getSelectedEpochTreeNodes{1};
recType = getRecordingTypeFromEpochList(parentNode.epochList);

highNode = parentNode.childBySplitValue('high').childBySplitValue(0).childBySplitValue(240);
minNode = parentNode.childBySplitValue('minimum').childBySplitValue(0).childBySplitValue(240);

temp_high = getMeanResponseTrace(highNode.epochList,recType);
temp_min = getMeanResponseTrace(minNode.epochList,recType);

figure; clf;
fig1=gca;
set(fig1,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Time (sec)')
set(get(fig1,'YLabel'),'String','V (mV)')

addLineToAxis(temp_high.timeVector,temp_high.mean,'high',fig1,'k','-','none')
addLineToAxis(temp_min.timeVector,temp_min.mean,'min',fig1,'r','-','none')
legend('9,000','Monitor off')


% % %
figure; clf;
fig2=gca;
set(fig2,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig2,'XLabel'),'String','Time (sec)')
set(get(fig2,'YLabel'),'String','V (mV)')
minNode = parentNode.childBySplitValue('minimum');
colors = pmkmp(minNode.children.length);
for cc = 1:minNode.children.length
    curNode = minNode.children(cc).childBySplitValue(480);
    temp = getMeanResponseTrace(curNode.epochList,recType);
    addLineToAxis(temp.timeVector,temp.mean,'high',fig2,colors(cc,:),'-','none')
    
end

