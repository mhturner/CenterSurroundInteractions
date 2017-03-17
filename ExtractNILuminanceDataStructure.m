%% CS NATURAL IMAGE LUMINANCE: tree

clear list; clear all; clear java; close all; clc; %#ok<CLJAVA,CLALL>
loader = edu.washington.rieke.Analysis.getEntityLoader();
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();
dataFolder = '/Users/maxturner/CurrentData/RFSurround/';
saveFileDirectory = '~/Documents/MATLAB/RFSurround/resources/SavedTreeFlags/';
import auimodel.*
import vuidocument.*
cd('~/Documents/MATLAB/RFSurround/')

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
    'protocolSettings(currentStimulus)'});
gui = epochTreeGUI(tree);

%% example tag at "exc" nodes. Select whole tree

NILuminanceData = struct;

parentNode = gui.getSelectedEpochTreeNodes{1};  
populationNodes = {};
ct = 0;
for nn = 1:parentNode.descendentsDepthFirst.length
    if parentNode.descendentsDepthFirst(nn).custom.get('isExample')
        ct = ct + 1;
        populationNodes(ct) = parentNode.descendentsDepthFirst(nn);
    end
end
for pp = 1:length(populationNodes)
    currentNode = populationNodes{pp};
    cellInfo = getCellInfoFromEpochList(currentNode.epochList);
    
    NILuminanceData(pp).cellID = cellInfo.cellID;
    NILuminanceData(pp).cellType = cellInfo.cellType;
    NILuminanceData(pp).NILuminance = [];
    for ii = 1:currentNode.childBySplitValue('CSNaturalImageLuminance').children.length
        pullNode = currentNode.childBySplitValue('CSNaturalImageLuminance').children(ii).childBySplitValue('Center');
        res = getMeanResponseTrace(pullNode.epochList,'exc');
        
        sampleRate = pullNode.epochList.firstValue.protocolSettings('sampleRate');
        fixationDuration = (pullNode.epochList.firstValue.protocolSettings('fixationDuration') / 1e3) *sampleRate; %data points
        startPoint = (pullNode.epochList.firstValue.protocolSettings('preTime') / 1e3) * sampleRate; %data points
        pad = ones(1,startPoint) .* pullNode.epochList.firstValue.protocolSettings('backgroundIntensity');
        stimTraceFix = convertJavaArrayList(pullNode.epochList.firstValue.protocolSettings('CenterIntensity'));
        stimTrace = [pad, kron(stimTraceFix,ones(1,fixationDuration)), pad];
        
        NILuminanceData(pp).NILuminance(ii).stimulus = stimTrace;
        NILuminanceData(pp).NILuminance(ii).response = res.mean;
    end
    noiseNode = currentNode.childBySplitValue('CenterSurroundNoise').children(1).childBySplitValue('Center');
    center = getLinearFilterAndPrediction(noiseNode.epochList,'exc',...
            'seedName','centerNoiseSeed');
    NILuminanceData(pp).Noise.stimulus = center.stimulus;
    NILuminanceData(pp).Noise.response = center.measuredResponse;
end

save(['NILuminanceData_',date],'NILuminanceData')
%%

figure(9); clf;
subplot(211);
plot(NILuminanceData(3).NILuminance(2).stimulus,'k')
subplot(212);
plot(NILuminanceData(3).NILuminance(2).response,'k')


