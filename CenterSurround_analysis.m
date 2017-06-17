%% RESET. NEW TREE.
clear list; clear all; clear java; close all; clc; %#ok<CLJAVA,CLALL>
loader = edu.washington.rieke.Analysis.getEntityLoader();
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();
dataFolder = '/Users/maxturner/CurrentData/RFSurround/';

saveFileDirectory = '~/Documents/MATLAB/RFSurround/resources/SavedTreeFlags/';
import auimodel.*
import vuidocument.*
cd('~/Documents/MATLAB/RFSurround/')

%% DOVES CS ADDITIVITY: tree
list = loader.loadEpochList([dataFolder,'DOVEScsAdditivity.mat'],dataFolder);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    'protocolSettings(stimulusIndex)',...
    recordingSplit_java,...
    'protocolSettings(currentStimulus)'});
gui = epochTreeGUI(tree);
%% DOVES CS ADDITIVITY: analysis and figs
% flag stimulusIndex nodes in population, example nodes are stimulusIndex
% select whole tree at root
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doDOVEScsAdditivityAnalysis(parentNode,...
    'exportFigs',true);


%% CENTER SURROUND NOISE: tree

list = loader.loadEpochList([dataFolder,'CenterSurroundNoise.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(useRandomSeed)',...
    'protocolSettings(currentStimulus)'});

%check that some parameters are consistent within recordings...
constantSettings = {'centerOffset','backgroundIntensity','preTime',...
    'tailTime','frameDwell','annulusInnerDiameter',...
    'centerDiameter','annulusOuterDiameter','noiseStdv'};
for cellTypeIndex = 1:tree.children.length
    for cellNameIndex = 1:tree.children(cellTypeIndex).children.length
        for recTypeIndex = 1:tree.children(cellTypeIndex).children(cellNameIndex).length
            newEL = tree.children(cellTypeIndex).children(cellNameIndex).children(recTypeIndex).epochList;
            [outEpochList, outParams] = makeUniformEpochList(newEL,constantSettings,[]);
        end
    end
end

gui = epochTreeGUI(tree);

%% CENTER SURROUND NOISE: cascade model plotting
% flag recType nodes in population, example nodes are useRandomSeed = 1
% select tree node (does ON & OFF in pop analysis together)
% eg Off is 20161129Ec1
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doCSLNAnalysis(parentNode,...
    'bins2D',15^2,... 
    'bins1D',20,...
    'exportFigs',true,...
    'convertToConductance',false,...
    'fitWithEquallyPopulatedBins',true);

%% CS NATURAL IMAGE LUMINANCE: tree
list = loader.loadEpochList([dataFolder,'CSNaturalImageLuminance.mat'],dataFolder);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

protocolSplit = @(list)splitOnShortProtocolID(list);
protocolSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, protocolSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    protocolSplit_java,...
    'protocolSettings(imageIndex)',...
    'protocolSettings(shuffleCenterSurround)',...
    'protocolSettings(currentStimulus)'});
gui = epochTreeGUI(tree);

%% CS NATURAL IMAGE LUMINANCE: analysis & figures
% flag recType nodes in population, example nodes are at imageIndex (under
%           CSNaturalImageLuminance node)
% select tree node (does ON & OFF in pop analysis together)
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};  
doCSNaturalImageLuminanceAnalysis(parentNode,...
    'exportFigs',true);


%% CORRELATED CS NOISE: tree
list = loader.loadEpochList([dataFolder,'CorrelatedCSNoise.mat'],dataFolder);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

protocolSplit = @(list)splitOnShortProtocolID(list);
protocolSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, protocolSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    protocolSplit_java,...
    'protocolSettings(csCorrelation)',...
    'protocolSettings(currentStimulus)'});
gui = epochTreeGUI(tree);

%% CORRELATED CS NOISE: analysis & figures
% flag recType nodes in population, example nodes are at CorrelatedCSNoise
%       node
% select tree node (does ON & OFF in pop analysis together)
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};  
doCorrelatedCSNoiseAnalysis(parentNode,...
    'exportFigs',true);

%% EXPANDING SPOTS: tree

list = loader.loadEpochList([dataFolder,'ExpandingSpots_ES.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    'protocolSettings(spotIntensity)',...
    recordingSplit_java,...
    'protocolSettings(currentSpotSize)'});
gui = epochTreeGUI(tree);
%% EXPANDING SPOTS: do analysis & make figs for exc vs. spikes analysis
%       Flag and example at spotIntensity. Select cell type
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doExcSpikesSurroundAnalysis(parentNode,'metric','integrated','exportFigs',false);
%% LINEAR EQUIVALENT DISC MOD SURROUND: tree
list = loader.loadEpochList([dataFolder,'LinearEquivalentDiscModSurround.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'currentPatchLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(imageName)',...
    locationSplit_java,...
    'protocolSettings(currentSurroundContrast)',...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% LINEAR EQUIVALENT DISC MOD SURROUND: analysis
% flag recType nodes, set patch location as examples
% select cell type as parentNode

clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doLEDModSurroundAnalysis(parentNode,...
    'metric','integrated','figureID','OFFspk');

% 
%% FLASHED GRATING MOD SURROUND: tree
list = loader.loadEpochList([dataFolder,'FlashedGratingModSurround.mat'],dataFolder);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

protocolSplit = @(list)splitOnShortProtocolID(list);
protocolSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, protocolSplit);


tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(gratingContrast)',...
    'protocolSettings(currentSurroundContrast)',...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% FLASHED GRATING MOD SURROUND: analysis & figs
% flag recType nodes, set grating contrast as examples
% select cell type as parentNode
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doFlashedGratingModSurroundAnalysis(parentNode,...
    'metric','integrated','exportFigs',true);


%% LINEAR EQUIVALENT DISC MIXED SURROUND: tree
list = loader.loadEpochList([dataFolder,'LEDMixedSurround.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'centerPatchLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

surroundSplit = @(list)splitOnJavaArrayList(list,'currentSurroundLocation');
surroundSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, surroundSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(imageName)',...
    locationSplit_java,...
    surroundSplit_java,...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% LINEAR EQUIVALENT DISC MIXED SURROUND: analysis
% flag recType nodes, set patch location as examples
% select cell type as parentNode

clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doLEDMixedSurroundAnalysis(parentNode,...
    'metric','integrated');

%% CONTRAST RESPONSE SPOTS: tree

list = loader.loadEpochList([dataFolder,'ContrastResponseSpots.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(maskDiameter)',...
    'protocolSettings(spotDiameter)',...
    'protocolSettings(currentSpotContrast)'});
gui = epochTreeGUI(tree);

%% CONTRAST RESPONSE SPOTS: flag population and set example(s)
%   Select cell type as parentNode
%       Flag spotDiameter node of cells in population    

% save flagged and example-d tree nodes for the future
saveFileID = 'CRF_HorizontalCell';
getFlaggedAndExampleNodes(tree, saveFileDirectory, saveFileID);
%% CONTRAST RESPONSE SPOTS: do analysis & make figs
%       Select cell type as parentNode
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doContrastResponseAnalysis(parentNode,'metric','peak',...
    'contrastPolarity',-1);


%% CS eye movement luminance
list = loader.loadEpochList([dataFolder,'CSEyeMovementLuminance.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

centerSplit = @(list)splitOnJavaArrayList(list,'centerOffset');
centerSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, centerSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,centerSplit_java,...
    'protocolSettings(stimulusIndex)',...
    'protocolSettings(currentStimulus)'});

gui = epochTreeGUI(tree);

%% CS eye movement luminance - simple plots
PSTHsigma = 20; %msec
parentNode = gui.getSelectedEpochTreeNodes{1};


centerResponse = getMeanResponseTrace(parentNode.childBySplitValue('Center').epochList,...
    'extracellular','PSTHsigma',PSTHsigma);
surroundResponse = getMeanResponseTrace(parentNode.childBySplitValue('Surround').epochList,...
    'extracellular','PSTHsigma',PSTHsigma);
csResponse = getMeanResponseTrace(parentNode.childBySplitValue('Center-Surround').epochList,...
    'extracellular','PSTHsigma',PSTHsigma);


figure(2); clf;
subplot(211); hold on;
plot(centerResponse.mean,'b')
plot(csResponse.mean,'k')
subplot(212); hold on;
plot(centerResponse.mean + surroundResponse.mean,'b')
plot(csResponse.mean,'k')


%% LINEAR EQUIVALENT CS ADDITIVITY: tree
list = loader.loadEpochList([dataFolder,'LinearEquivalentCSAdditivity.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'currentPatchLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(imageName)',...
    'protocolSettings(currentCenter)',...
    'protocolSettings(currentSurround)',...
    locationSplit_java});

gui = epochTreeGUI(tree);

%% flag recType nodes and set an example ON and OFF image node
% select whole tree
clc;
parentNode = gui.getSelectedEpochTreeNodes{1};
doLECSAnalysis(parentNode,'metric','integrated',...
    'figureID','OFFparSpikes');

%%
clc;
load('NaturalImageFlashLibrary_101716.mat')
imageNames = fieldnames(imageData);
resourcesDir = '~/Documents/MATLAB/turner-package/resources';
patchSize = 100; %pixels
MicronsPerPixel = 6.6;
patchesPerImage = 100;
ImageID = 'imk01769';

RFmodel = ThreeLayerReceptiveFieldModel;
RFmodel.MicronsPerPixel = MicronsPerPixel;
RFmodel.makeRfComponents(patchSize);

ImageX = 1536; ImageY = 1024;

fileId=fopen([resourcesDir, '/VHsubsample_20160105', '/', ImageID,'.iml'],'rb','ieee-be');
img = fread(fileId, [1536,1024], 'uint16');
contrastImg = img - mean(img(:));
contrastImg = contrastImg ./ max(contrastImg(:));

noBins = 2.*patchesPerImage; %from no. image patches to show
[N, edges, bin] = histcounts(imageData.(ImageID).responseDifferences,noBins);
populatedBins = unique(bin);
%pluck one patch from each bin
pullInds = arrayfun(@(b) find(b == bin,1),populatedBins);
location = imageData.(ImageID).location(pullInds,:);

tempResponse = struct;
tempResponse.SubunitCenterPlusLinearSurround = [];
tempResponse.SubunitCenterPlusSubunitSurround = [];
for rr = 1:patchesPerImage
    x = location(rr,1); y = location(rr,2);
    % get patch
    ContrastPatch = -contrastImg(x-patchSize/2+1:x+patchSize/2,y-patchSize/2+1:y+patchSize/2);
    % get RF model responses
    responseStructure = RFmodel.getResponse(ContrastPatch);
    tempResponse.SubunitCenterPlusLinearSurround(rr) = responseStructure.CenterSurround.SharedNonlinearity;
    tempResponse.SubunitCenterPlusSubunitSurround(rr) = responseStructure.CenterSurround.LNLNSamePolarity;
end

figure(3); clf;
fig2=gca;
set(fig2,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig2,'XLabel'),'String','Center + nonlinear surround')
set(get(fig2,'YLabel'),'String','Center + linear surround')

colors = pmkmp(3);


addLineToAxis(tempResponse.SubunitCenterPlusSubunitSurround,tempResponse.SubunitCenterPlusLinearSurround,...
    'modelResp',fig2,colors(2,:),'none','o')
addLineToAxis([0 0.3],[0 0.3],...
    'unity',fig2,'k',':','none')

% makeAxisStruct(fig2,'LECS_model' ,'COSYNE2017Figs')

%% annulus stim images for figures...
load('NaturalImageFlashLibrary_101716.mat')
imageNames = fieldnames(imageData);
resourcesDir = '~/Documents/MATLAB/turner-package/resources';
ImageID = 'imk00152';
fileId=fopen([resourcesDir, '/VHsubsample_20160105', '/', ImageID,'.iml'],'rb','ieee-be');
img = fread(fileId, [1536,1024], 'uint16');
img = 255 .* img ./ max(img(:));
patchSize = 200; %pixels
MicronsPerPixel = 3.3;

patchesPerImage = 10;

rr = 9;

noBins = 2.*patchesPerImage; %from no. image patches to show
[N, edges, bin] = histcounts(imageData.(ImageID).responseDifferences,noBins);
populatedBins = unique(bin);
%pluck one patch from each bin
pullInds = arrayfun(@(b) find(b == bin,1),populatedBins);
location = imageData.(ImageID).location(pullInds,:);
    
x = location(rr,1); y = location(rr,2);

x = 1284; y = 657;

% get patch
patch = img(x-patchSize/2+1:x+patchSize/2,y-patchSize/2+1:y+patchSize/2);
origPatch = patch;

centerRadius = round(100/MicronsPerPixel);
surroundInner = round(150/MicronsPerPixel);
surroundOuter = round(300/MicronsPerPixel);

[rr, cc] = meshgrid(1:patchSize,1:patchSize);
tempDist = sqrt((rr-(patchSize/2)).^2+(cc-(patchSize/2)).^2);
centerBinary = ~(tempDist < centerRadius);
surroundBinary = ~(tempDist < surroundOuter & tempDist > surroundInner);
binary = min(centerBinary,surroundBinary);

patch(binary) = mean(img(:));
patchCS = patch;
patchCL = patch;
patchCL(~surroundBinary) = mean(patchCS(~surroundBinary));

patchC = origPatch;
patchS = origPatch;
patchC(centerBinary) = mean(img(:));
patchS(surroundBinary) = mean(img(:));

figure(2); clf;
subplot(211)
image(patchCS'); colormap(gray); axis image; axis off;
caxis([min(img(:)) max(img(:))])
subplot(212)
image(patchCL'); colormap(gray); axis image; axis off;
caxis([min(img(:)) max(img(:))])

figure(3); clf;
subplot(311)
image(patchC'); colormap(gray); axis image; axis off;
caxis([min(img(:)) max(img(:))])
subplot(312)
image(patchS'); colormap(gray); axis image; axis off;
caxis([min(img(:)) max(img(:))])
subplot(313)
image(patchCS'); colormap(gray); axis image; axis off;
caxis([min(img(:)) max(img(:))])

%% % % % % Linearity of center modulated by surround % % % % % % 
colors = hsv(length(allResp) + 1);
for im = 1:length(allResp);
    color = colors(im,:);
    responseMatrix = allResp{im};
figure(2);
subplot(131); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,1), responseMatrix(:,4),'Marker','o','Color',color,'LineStyle','none')
xlabel('Image'); ylabel('Disc')


subplot(132); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,7), responseMatrix(:,6),'Marker','o','Color',color,'LineStyle','none')
xlabel('Image + equiv surround'); ylabel('Disc + equiv surround')
title('Linearity of center modulated by surround')

subplot(133); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,3), responseMatrix(:,8),'Marker','o','Color',color,'LineStyle','none')
xlabel('Image + image surround'); ylabel('Disc + image surround')


% subplot(235); plot([0 30],[0 30],'k--'); hold on;
% plot(responseMatrix(:,3), responseMatrix(:,6),'Marker','o','Color',color,'LineStyle','none')
% xlabel('Image + image surround'); ylabel('Disc + equiv surround')
% title('Matched')
% 
% subplot(236); plot([0 30],[0 30],'k--'); hold on;
% plot(responseMatrix(:,7), responseMatrix(:,8),'Marker','o','Color',color,'LineStyle','none')
% xlabel('Image + equiv surround'); ylabel('Disc + image surround')
% title('Anti-matched')

end

%% % % % % % Linearity of surround % % % % % % % % % 
for im = 1:length(allResp);
color = colors(im,:);
responseMatrix = allResp{im};
figure(3);
% subplot(131); plot([0 30],[0 30],'k--'); hold on;
% plot(responseMatrix(:,2), responseMatrix(:,5),'Marker','o','Color',color,'LineStyle','none')
% xlabel('image surround'); ylabel('equiv surround')


subplot(121); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,8), responseMatrix(:,6),'Marker','o','Color',color,'LineStyle','none')
xlabel('Disc + image surround'); ylabel('Disc + equiv surround')
title('Linearity of surround')

subplot(122); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,3), responseMatrix(:,7),'Marker','o','Color',color,'LineStyle','none')
xlabel('Image + image surround'); ylabel('Image + equiv surround')
end
%% % % % % Center-Surround Additivity % % % % % % % % % 
for im = 1:length(allResp);
color = colors(im,:);
responseMatrix = allResp{im};
figure(4);
subplot(221); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,1)+responseMatrix(:,2), responseMatrix(:,3),'Marker','o','Color',color,'LineStyle','none')
xlabel('R(C) + R(S)'); ylabel('R(C + S)')
title('Image center, Image surround')

subplot(222); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,4)+responseMatrix(:,5), responseMatrix(:,6),'Marker','o','Color',color,'LineStyle','none')
xlabel('R(C) + R(S)'); ylabel('R(C + S)')
title('Equiv center, Equiv surround')

subplot(223); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,4)+responseMatrix(:,2), responseMatrix(:,8),'Marker','o','Color',color,'LineStyle','none')
xlabel('R(C) + R(S)'); ylabel('R(C + S)')
title('Equiv center, Image surround')

subplot(224); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,1)+responseMatrix(:,5), responseMatrix(:,7),'Marker','o','Color',color,'LineStyle','none')
xlabel('R(C) + R(S)'); ylabel('R(C + S)')
title('Image center, Equiv surround')
end

%% % % % % Gain control by surround % % % % % % % % % 
for im = 1:length(allResp);
color = colors(im,:);
responseMatrix = allResp{im};
figure(6);
subplot(221); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,1), responseMatrix(:,3),'Marker','o','Color',color,'LineStyle','none')
xlabel('C'); ylabel('C + S')
title('Image center, Image surround')

subplot(222); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,4), responseMatrix(:,6),'Marker','o','Color',color,'LineStyle','none')
xlabel('C'); ylabel('C  + S')
title('Equiv center, Equiv surround')

subplot(223); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,1), responseMatrix(:,7),'Marker','o','Color',color,'LineStyle','none')
xlabel('C'); ylabel('C + S')
title('Image center, Equiv surround')

subplot(224); plot([0 30],[0 30],'k--'); hold on;
plot(responseMatrix(:,4), responseMatrix(:,8),'Marker','o','Color',color,'LineStyle','none')
xlabel('C'); ylabel('C + S')
title('Equiv center, Image surround')
end

%% TEXTURE CS ADDITIVITY: tree
list = loader.loadEpochList([dataFolder,'TextureCSAdditivity.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(contrast)',...
    'protocolSettings(currentStimulus)'});

gui = epochTreeGUI(tree);

%% select contrast node
 
parentNode = gui.getSelectedEpochTreeNodes{1};

treeC = riekesuite.analysis.buildTree(parentNode.childBySplitValue('Center').epochList,...
    {'protocolSettings(currentCenterSigma)'});
treeS = riekesuite.analysis.buildTree(parentNode.childBySplitValue('Surround').epochList, ...
    {'protocolSettings(currentSurroundSigma)'});
treeCS = riekesuite.analysis.buildTree(parentNode.childBySplitValue('Center-Surround').epochList, ...
    {'protocolSettings(currentCenterSigma)','protocolSettings(currentSurroundSigma)'});


responses.Center = [];
sigmas.Center = [];
for cc = 1:treeC.children.length
    sigmas.Center(cc) = treeC.children(cc).splitValue;
    newEpochList = treeC.children(cc).epochList;
    newResp = getResponseAmplitudeStats(newEpochList,'extracellular');
    responses.Center(cc) = newResp.integrated.mean;
end

responses.Surround = [];
sigmas.Surround = [];
for ss = 1:treeS.children.length
    sigmas.Surround(ss) = treeS.children(ss).splitValue;
    newEpochList = treeS.children(ss).epochList;
    newResp = getResponseAmplitudeStats(newEpochList,'extracellular');
    responses.Surround(ss) = newResp.integrated.mean;
end

responses.CenterSurround = [];
for cc = 1:treeCS.children.length
    for ss = 1:treeCS.children(cc).children.length
        newEpochList = treeCS.children(cc).children(ss).epochList;
        newResp = getResponseAmplitudeStats(newEpochList,'extracellular');
        responses.CenterSurround(ss,cc) = newResp.integrated.mean;
    end
end

%%
colors = pmkmp(3);
figure(2); clf;
subplot(211); hold on
plot(sigmas.Center,responses.Center,'Color',colors(1,:),'Marker','o');
plot(sigmas.Surround,responses.Surround,'Color',colors(2,:),'Marker','o');
legend('Center','Surround')
xlabel('Spatial scale sigma (um)')
ylabel('Response (spikes)')
subplot(212)
pcolor(sigmas.Center,sigmas.Surround,responses.CenterSurround)
colormap(hot)
xlabel('Center'); ylabel('Surround'); colorbar;

%% CenterF2 Plus Surround: tree

list = loader.loadEpochList([dataFolder,'CenterF2PlusSurround.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(flashDelay)',...
    'protocolSettings(flashDuration)',...
    'protocolSettings(temporalFrequency)',...
    'protocolSettings(centerContrast)',...
    'protocolSettings(surroundContrast)',...
    'protocolSettings(currentStimulus)'});
gui = epochTreeGUI(tree);

%%

parentNode = gui.getSelectedEpochTreeNodes{1};

surroundAloneNode = parentNode.childBySplitValue(0).childBySplitValue(0.9);
surroundMean = getMeanResponseTrace(surroundAloneNode.epochList,'exc');

TestNode = parentNode.childBySplitValue(0.25).childBySplitValue(0.9);
F1mean = getMeanResponseTrace(TestNode.childBySplitValue('FullField').epochList,'exc');
F2mean = getMeanResponseTrace(TestNode.childBySplitValue('SplitField').epochList,'exc');


figure(2); clf;
subplot(211); hold on; 
plot(surroundMean.timeVector,surroundMean.mean,'k')
plot(F1mean.timeVector,F1mean.mean,'b')
plot(F2mean.timeVector,F2mean.mean,'r')

subplot(212); hold on;
plot(surroundMean.timeVector,F1mean.mean-surroundMean.mean,'b')
plot(surroundMean.timeVector,F2mean.mean-surroundMean.mean,'r')

periodLen = round(1 / TestNode.epochList.firstValue.protocolSettings('temporalFrequency') .* 1e4); %dataPts
prePoints = TestNode.epochList.firstValue.protocolSettings('preTime') .* 10; %dataPts
startPoint = prePoints + 2*periodLen;
endPoint = startPoint + 5.*periodLen;
cycles = zeros(5,periodLen);
cyclesF2 = zeros(5,periodLen);
for cc = 1:5
    currentStart = startPoint + (cc-1)*periodLen + 1;
    currentEnd = startPoint + cc * periodLen;
    cycles(cc,:) = F1mean.mean(currentStart:currentEnd);
    cyclesF2(cc,:) = F2mean.mean(currentStart:currentEnd);
end

cycleAvg = -mean(cycles,1); %flip exc to positive
F2Resp = -mean(cyclesF2,2);
timeVec = (0:(periodLen-1)) ./ 1e4; %sec
tF = 8;

% b(1) is amplitude. Freq is fixed. Offset & phase allowed to
% vary
modelFun = @(b,t)(b(1).*(sin(2*pi*t.*tF + b(2))) + b(3));
beta0 = [max(cycleAvg);0;0];
beta = nlinfit(timeVec,cycleAvg,modelFun,beta0);

stimulus_F1 = modelFun([0.5,beta(2),0],timeVec);
stimulus_RightF2 = modelFun([0.5,beta(2),0],timeVec);
stimulus_LeftF2 = modelFun([0.5,beta(2)+pi,0],timeVec);

params0 = [1e4, 0.5, -4, 0];
fitRes = fitCRF_cumGauss(stimulus_F1,cycleAvg,params0);

figure(4); clf;
subplot(121)
plot(timeVec,cycleAvg,'b')
subplot(122); hold on;
plot(stimulus_F1,cycleAvg,'b')
fitXX = -1:0.01:1;
fitYY = CRFcumGauss(fitXX,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);
plot(fitXX, fitYY,'k')

%long trace predictions
timeVec = (0:(16*periodLen-1)) ./ 1e4; %sec
stimulus_F1 = modelFun([0.9,beta(2),0],timeVec);
stimulus_RightF2 = modelFun([0.9,beta(2),0],timeVec);
stimulus_LeftF2 = modelFun([0.9,beta(2)+pi,0],timeVec);

stepTime = find(timeVec >= 0.8625,1);
stepStimulus = zeros(size(timeVec));
stepStimulus(stepTime + 1 : stepTime + 2500) = 0.5 * 0.9;

%indep stims...
F1prediction = CRFcumGauss(stimulus_F1,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

F2left = 0.5*CRFcumGauss(stimulus_LeftF2,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

F2right = 0.5*CRFcumGauss(stimulus_RightF2,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

stepPrediction = CRFcumGauss(stepStimulus,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

F2prediction = F2left + F2right;

figure(5); clf;
hold on;
plot(timeVec,F1prediction,'b')
plot(timeVec,F2prediction,'r')
plot(timeVec,stepPrediction,'k')


%plus step stims
F1prediction = CRFcumGauss(stimulus_F1 + stepStimulus,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

F2left = 0.5*CRFcumGauss(stimulus_LeftF2 + stepStimulus,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

F2right = 0.5*CRFcumGauss(stimulus_RightF2 + stepStimulus,...
    fitRes.alphaScale,...
    fitRes.betaSens,...
    fitRes.gammaXoff,...
    fitRes.epsilonYoff);

F2prediction = F2left + F2right;

figure(6); clf;
hold on;
plot(timeVec,F1prediction,'b')
plot(timeVec,F2prediction,'r')
plot(timeVec,stepPrediction,'k')



