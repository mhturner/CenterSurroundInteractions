%% RESET. NEW TREE.
clear list; clear all; clear java; close all; clc; %#ok<CLJAVA,CLALL>
loader = edu.washington.rieke.Analysis.getEntityLoader();
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();
dataFolder = '/Users/mhturner/Dropbox/CurrentData/RFSurround/';

saveFileDirectory = '~/Dropbox/RiekeLab/Analysis/MATLAB/RFSurround/resources/SavedTreeFlags/';
import auimodel.*
import vuidocument.*
cd('~/Dropbox/RiekeLab/Analysis/MATLAB/RFSurround/')

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
% eg OFF 20170411Ec4_13
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
% eg On is 20170822Ec2
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doCSLNAnalysis(parentNode,...
    'bins2D',15^2,... 
    'bins1D',20,...
    'exportFigs',false,...
    'convertToConductance',true,...
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
% Eg off 20170214Ec5 (im 5)
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
% eg cell 
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
% Eg off spike: 20170110Ec2

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

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(gratingContrast)',...
    'protocolSettings(currentSurroundContrast)',...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% FLASHED GRATING MOD SURROUND: analysis & figs
% flag recType nodes, set grating contrast as examples
% select cell type as parentNode
% OFF eg = 20170214Ec5, 0.9 contrast
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doFlashedGratingModSurroundAnalysis(parentNode,...
    'metric','integrated','exportFigs',true);

%% FLASHED GRATING CORRELATED SURROUND: tree
list = loader.loadEpochList([dataFolder,'FlashedGratingCorrelatedSurround.mat'],dataFolder);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(gratingContrast)',...
    'protocolSettings(currentIntensity)',...
    'protocolSettings(surroundTag)',...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% FLASHED GRATING CORRELATED SURROUND: analysis & figs
% flag recType nodes, set rec type as example cells
% select cell type as parentNode
% eg Off: 20170925Ec3
clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doFlashedGratingCorrelatedSurroundAnalysis(parentNode,...
    'metric','integrated','exportFigs',true);


%% LINEAR EQUIVALENT DISC MIXED SURROUND: tree
list = loader.loadEpochList([dataFolder,'LEDMixedSurround.mat'],dataFolder);

recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

locationSplit = @(list)splitOnJavaArrayList(list,'currentCenterLocation');
locationSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, locationSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(imageName)',...
    locationSplit_java,...
    'protocolSettings(surroundTag)',...
    'protocolSettings(stimulusTag)'});

gui = epochTreeGUI(tree);

%% LINEAR EQUIVALENT DISC MIXED SURROUND: analysis
% flag recType nodes
% set image ID as example. E.g. = 
% select cell type as parentNode

clc; CloseAllFiguresExceptGUI();
parentNode = gui.getSelectedEpochTreeNodes{1};
doLEDMixedSurroundAnalysis(parentNode,...
    'metric','integrated','exportFigs',false);

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





