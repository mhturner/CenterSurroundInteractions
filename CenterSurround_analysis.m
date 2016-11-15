%% RESET. NEW TREE.
clear list; clear all; clear java; close all; clc; %#ok<CLJAVA,CLALL>
loader = edu.washington.rieke.Analysis.getEntityLoader();
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();
dataFolder = '/Users/maxturner/CurrentData/RFSurround/';
saveFileDirectory = '~/Documents/MATLAB/RFSurround/resources/SavedTreeFlags/';
import auimodel.*
import vuidocument.*
cd('~/Documents/MATLAB/RFSurround/')

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

%% CENTER SURROUND NOISE: tree


list = loader.loadEpochList([dataFolder,'CenterSurroundNoise.mat'],dataFolder);
recordingSplit = @(list)splitOnRecKeyword(list);
recordingSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, recordingSplit);

cellTypeSplit = @(list)splitOnCellType(list);
cellTypeSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, cellTypeSplit);

tree = riekesuite.analysis.buildTree(list, {cellTypeSplit_java,'cell.label',...
    recordingSplit_java,...
    'protocolSettings(currentStimulus)'});

%check that some parameters are consistent within recordings...
constantSettings = {'centerOffset','backgroundIntensity','preTime',...
    'stimTime','tailTime','frameDwell','annulusInnerDiameter',...
    'centerDiameter','annulusOuterDiameter','useRandomSeed','noiseStdv'};
for cellTypeIndex = 1:tree.children.length
    for cellNameIndex = 1:tree.children(cellTypeIndex).children.length
        for recTypeIndex = 1:tree.children(cellTypeIndex).children(cellNameIndex).length
            newEL = tree.children(cellTypeIndex).children(cellNameIndex).children(recTypeIndex).epochList;
            [outEpochList, outParams] = makeUniformEpochList(newEL,constantSettings,[]);
        end
    end
end

gui = epochTreeGUI(tree);

%% CENTER SURROUND NOISE: example cell plotting
% choose rec type node
clc

numberOfBins = 8^2;

parentNode = gui.getSelectedEpochTreeNodes{1};
cellInfo = getCellInfoFromEpochList(parentNode.epochList);
recType = getRecordingTypeFromEpochList(parentNode.epochList);

% % % % % % % % GET LINEAR FILTERS AND GENERATOR SIGNALS % % % % % % % % % % % % % % % % 
% center:
center = getLinearFilterAndPrediction(parentNode.childBySplitValue('Center').epochList,recType,...
    'seedName','centerNoiseSeed','numberOfBins',numberOfBins);
% surround:
surround = getLinearFilterAndPrediction(parentNode.childBySplitValue('Surround').epochList,recType,...
    'seedName','surroundNoiseSeed','numberOfBins',numberOfBins);
% center + surround: just to test against measured response
centerSurround = getLinearFilterAndPrediction(parentNode.childBySplitValue('Center-Surround').epochList,recType,...
    'seedName','centerNoiseSeed','numberOfBins',numberOfBins);

figure(2); clf;
subplot(321)
plot(center.filterTimeVector,center.LinearFilter,'b');
xlabel('Time (s)'); title('Center')
subplot(322)
plot(surround.filterTimeVector,surround.LinearFilter,'r');
xlabel('Time (s)'); title('Surround')

% % % % % % % % MODEL FITTING % % % % % % % % % % % % % % % % 
epochLen = length(centerSurround.measuredResponse) / centerSurround.n;
epochToHold = 5; %epoch number to hold out
testDataInds = ((epochToHold-1)*epochLen + 1):(epochToHold*epochLen);
fitDataInds = setdiff(1:length(centerSurround.measuredResponse),testDataInds);

%bin up and shape training data:
[~,centerGS,~,centerBinID] = ...
    histcounts_equallyPopulatedBins(center.generatorSignal(fitDataInds),sqrt(numberOfBins));
[~,surroundGS,~,surroundBinID] = ...
    histcounts_equallyPopulatedBins(surround.generatorSignal(fitDataInds),sqrt(numberOfBins));
responseMean = zeros(sqrt(numberOfBins));
for xx = 1:sqrt(numberOfBins)
    for yy = 1:sqrt(numberOfBins)
        jointInds = intersect(find(centerBinID == xx),find(surroundBinID == yy));
        tempResp = centerSurround.measuredResponse(fitDataInds);
        responseMean(yy,xx) = mean(tempResp(jointInds));
    end
end
%Fit Joint nonlinearity model:
% params = [alpha mu1 mu2 std1 std2 corr12 epsilon]
params0 = [max(responseMean(:)), 400, 0,50, 50,0 0];
fitRes_joint = fitNLinearity_2D(centerGS,surroundGS,responseMean,params0);

cc = linspace(min(centerGS),max(centerGS),20);
ss = linspace(min(surroundGS),max(surroundGS),20);
[CC,SS] = meshgrid(cc',ss');
fitSurface = JointNLin_mvcn(CC(:)',SS(:)',fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);
fitSurface = reshape(fitSurface,length(ss),length(cc));

figure(2);
subplot(325); hold on;
stem3(centerGS,surroundGS,responseMean,'filled');
xlabel('Center'); ylabel('Surround'); zlabel('Response (pA)')
meshc(CC,SS,fitSurface)
title('Joint')
view(4,10);

% Fit Independent nonlinearity model:
%params is [alphaC, betaC, gammaC,...
%           alphaS, betaS, gammaS, epsilon]
params0 = [max(responseMean(:)), 0.002, -0.5,...
    max(responseMean(:)), 0.01, -1,...
    -400];
fitRes_indep = fitCSModel_IndependentNL(centerGS,surroundGS,responseMean,params0);

fitSurface = CSModel_IndependentNL(CC(:)',SS(:)',...
    fitRes_indep.alphaC,fitRes_indep.betaC,...
    fitRes_indep.gammaC,fitRes_indep.alphaS,...
    fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);

fitSurface = reshape(fitSurface,length(ss),length(cc));

figure(2);
subplot(323); hold on;
plot(center.nonlinearity.fitXX,center.nonlinearity.fitYY,'b')
plot(surround.nonlinearity.fitXX,surround.nonlinearity.fitYY,'r')

figure(2);
subplot(324); hold on;
stem3(centerGS,surroundGS,responseMean,'filled');
xlabel('Center'); ylabel('Surround'); zlabel('Response (pA)')
meshc(CC,SS,fitSurface)
title('Indep.')
view(4,10);

% Fit Shared nonlinearity model:
% params is [a, alpha, beta, gamma, epsilon]
params0=[2, max(responseMean(:)), 0.01, 0, 0]';
fitRes_shared = fitCSModel_SharedNL(centerGS,surroundGS,responseMean,params0);

fitSurface = CSModel_SharedNL(CC(:)',SS(:)',...
    fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
    fitRes_shared.gamma,fitRes_shared.epsilon);
fitSurface = reshape(fitSurface,length(ss),length(cc));

figure(2);
subplot(326); hold on;
stem3(centerGS,surroundGS,responseMean,'filled');
xlabel('Center'); ylabel('Surround'); zlabel('Response (pA)')
meshc(CC,SS,fitSurface)
title('Shared')
view(4,10);

% % % % % % % % PREDICTIONS % % % % % % % % % % % % % % % % 
centerGS = center.generatorSignal(testDataInds);
surroundGS = surround.generatorSignal(testDataInds);
measuredResponse = centerSurround.measuredResponse(testDataInds);

centerAloneResponse = center.measuredResponse(testDataInds);
surroundAloneResponse = surround.measuredResponse(testDataInds);

linearSummedResponse = centerAloneResponse + surroundAloneResponse;

predictedResponse_joint = JointNLin_mvcn(centerGS, surroundGS,...
    fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);

predictedResponse_indep = CSModel_IndependentNL(centerGS, surroundGS,...
    fitRes_indep.alphaC,fitRes_indep.betaC,...
    fitRes_indep.gammaC,fitRes_indep.alphaS,...
    fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);

predictedResponse_shared = CSModel_SharedNL(centerGS, surroundGS,...
    fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
    fitRes_shared.gamma,fitRes_shared.epsilon);

%Model predictions:
figure(3); clf;
subplot(331); hold on;
plot(predictedResponse_joint, measuredResponse, 'go')
cc = corr(predictedResponse_joint, measuredResponse');
title(num2str(cc))
xlabel('Predicted'); ylabel('Measured');
subplot(3,3,2:3); hold on
plot(predictedResponse_joint, 'g')
plot(measuredResponse, 'k')
title('Joint')

subplot(334); hold on;
plot(predictedResponse_indep, measuredResponse, 'bo')
cc = corr(predictedResponse_indep', measuredResponse');
title(num2str(cc))
xlabel('Predicted'); ylabel('Measured');
subplot(3,3,5:6); hold on
plot(predictedResponse_indep, 'b')
plot(measuredResponse, 'k')
title('Indep')

subplot(337); hold on;
plot(predictedResponse_shared, measuredResponse, 'ro')
cc = corr(predictedResponse_shared', measuredResponse');
title(num2str(cc))
xlabel('Predicted'); ylabel('Measured');
subplot(3,3,8:9); hold on
plot(predictedResponse_shared, 'r')
plot(measuredResponse, 'k')
title('Shared')

%CS versus response to c or s alone:
figure(4); clf;
subplot(331); hold on;
plot(centerAloneResponse, measuredResponse, 'go')
cc = corr(centerAloneResponse', measuredResponse');
title(num2str(cc))
xlabel('R(C)'); ylabel('R(C + S)');
subplot(3,3,2:3); hold on
plot(centerAloneResponse, 'g')
plot(measuredResponse, 'k')
title('Center alone')

subplot(334); hold on;
plot(surroundAloneResponse, measuredResponse, 'bo')
cc = corr(surroundAloneResponse', measuredResponse');
title(num2str(cc))
xlabel('R(S)'); ylabel('R(C + S)');
subplot(3,3,5:6); hold on
plot(surroundAloneResponse, 'b')
plot(measuredResponse, 'k')
title('Surround alone')

subplot(337); hold on;
plot(linearSummedResponse, measuredResponse, 'ro')
cc = corr(linearSummedResponse', measuredResponse');
title(num2str(cc))
xlabel('R(C) + R(S)'); ylabel('R(C + S)');
subplot(3,3,8:9); hold on
plot(linearSummedResponse, 'r')
plot(measuredResponse, 'k')
title('Linear sum')

%% Slices through joint nonlinearity: does center nlin change with surround?
[~,centerGS,~,centerBinID] = ...
    histcounts_equallyPopulatedBins(center.generatorSignal(fitDataInds),sqrt(numberOfBins));
[~,surroundGS,~,surroundBinID] = ...
    histcounts_equallyPopulatedBins(surround.generatorSignal(fitDataInds),sqrt(numberOfBins));
responseMean = zeros(sqrt(numberOfBins));
for xx = 1:sqrt(numberOfBins)
    for yy = 1:sqrt(numberOfBins)
        jointInds = intersect(find(centerBinID == xx),find(surroundBinID == yy));
        tempResp = centerSurround.measuredResponse(fitDataInds);
        responseMean(yy,xx) = mean(tempResp(jointInds));
    end
end


figure(7); clf;
subplot(211); hold on;
p.z = plot([0 0.5],[0 0],'k--');
p.c = plot(center.filterTimeVector,center.LinearFilter,'b');
xlabel('Time (s)'); xlim([0 0.5]); title('Center');
subplot(212); hold on;
p.z = plot([0 0.5],[0 0],'k--');
p.s = plot(surround.filterTimeVector,surround.LinearFilter,'r');
xlabel('Time (s)'); xlim([0 0.5]); title('Surround');

figure(8); clf; subplot(211)
mesh(centerGS,surroundGS,responseMean);
xlabel('Center'); ylabel('Surround'); zlabel('Response (pA)')

figure(8); subplot(212); hold on;
plot([0 0],[min(responseMean(:)) max(responseMean(:))],'k--')
plot([min(centerGS(:)) max(centerGS(:))],[0 0], 'k--')
colors = jet(length(surroundGS));
for ss = 1:length(surroundGS)
    p.(['s',num2str(ss)]) = plot(centerGS,responseMean(ss,:),'Color',colors(ss,:));
end

xlim([min(centerGS(:)) max(centerGS(:))]);
ylim([min(responseMean(:)) max(responseMean(:))]);

xlabel('Center activation'); ylabel('Response (pA)')
legend([p.s1, p.s8], ['Surround ',num2str(surroundGS(1))],['Surround ',num2str(surroundGS(8))])
%%

figure(9); clf; subplot(211); hold on;
colors = jet(length(surroundGS));
plot([0 0],[min(responseMean(:)) max(responseMean(:))],'k--')
plot([min(centerGS + surroundGS) max(centerGS + surroundGS)],[0 0], 'k--')
for ii = 1:8
   currentX = centerGS(ii) + surroundGS;
   currentY = responseMean(:,ii);
    plot(currentX,currentY,'Color',colors(ii,:));
end
ylabel('Response'); xlabel('Center + Modulated Surround')
xlim([min(centerGS + surroundGS) max(centerGS + surroundGS)]);
ylim([min(responseMean(:)) max(responseMean(:))]);

subplot(212); hold on;
colors = jet(length(surroundGS));
plot([0 0],[min(responseMean(:)) max(responseMean(:))],'k--')
plot([min(centerGS + surroundGS) max(centerGS + surroundGS)],[0 0], 'k--')
for ii = 1:8
   currentX = surroundGS(ii) + centerGS;
   currentY = responseMean(ii,:);
    plot(currentX,currentY,'Color',colors(ii,:));
end
ylabel('Response'); xlabel('Surround + Modulated Center')
xlim([min(centerGS + surroundGS) max(centerGS + surroundGS)]);
ylim([min(responseMean(:)) max(responseMean(:))]);
%% all data: model-free, where are failures of additivity?
numberOfBins = 8^2;

measuredResponse = centerSurround.measuredResponse;

centerAloneResponse = center.measuredResponse;
surroundAloneResponse = surround.measuredResponse;
linearSummedResponse = centerAloneResponse + surroundAloneResponse;

err = measuredResponse - linearSummedResponse;

%bin up and shape data:
[~,centerGS,~,centerBinID] = ...
    histcounts_equallyPopulatedBins(center.generatorSignal,sqrt(numberOfBins));
[~,surroundGS,~,surroundBinID] = ...
    histcounts_equallyPopulatedBins(surround.generatorSignal,sqrt(numberOfBins));
errorMatrix = zeros(sqrt(numberOfBins));
for xx = 1:sqrt(numberOfBins)
    for yy = 1:sqrt(numberOfBins)
        jointInds = intersect(find(centerBinID == xx),find(surroundBinID == yy));
        errorMatrix(yy,xx) = mean(err(jointInds));
    end
end

figure(5); clf;
subplot(321)
plot(center.filterTimeVector,center.LinearFilter,'b')
hold on; plot([0 0.8],[0 0],'k:')
xlabel('Time (s)'); title('Center')
xlim([0 0.6]); ylim([min(center.LinearFilter) max(center.LinearFilter)])

subplot(322)
plot(surround.filterTimeVector,surround.LinearFilter,'r')
hold on; plot([0 0.8],[0 0],'k:')
xlabel('Time (s)'); title('Surround')
xlim([0 0.6]); ylim([min(surround.LinearFilter) max(surround.LinearFilter)])


subplot(3,2,3:6)
surfc(centerGS,surroundGS,errorMatrix);
% axis square; view(-16,22)
axis square; view(164,32)
xlabel('Center gen. signal'); ylabel('Surround gen. signal');
zlabel('Measured - Linear sum');

%% eye movement luminace trajectories: histogram in c/s gen signal space
load('~/Documents/MATLAB/Analysis/NatImages/Stimuli/SaccadeLuminanceTrajectoryStimuli_20160919.mat')

numberOfBins_em = 100^2;
 
centerGenSignal = [];
surroundGenSignal = [];
filterSampleFreq = length(center.filterTimeVector) / center.filterTimeVector(end); %Hz
for ss = 1:length(luminanceData)
    cStim = resample(luminanceData(ss).centerTrajectory,30,200  );
    cStim = (cStim) ./ luminanceData(ss).ImageMax; %stim as presented

    %convert to contrast (relative to mean) for filter convolution
    imMean = (luminanceData(ss).ImageMean  ./ luminanceData(ss).ImageMax);
    cStim = (cStim - imMean) / imMean;
    
    linearPrediction = conv(cStim,center.LinearFilter);
    linearPrediction = linearPrediction(1:length(cStim));
    centerGenSignal = cat(2,centerGenSignal,linearPrediction);

    sStim = resample(luminanceData(ss).surroundTrajectory,30,200);
    
    sStim = (sStim) ./ luminanceData(ss).ImageMax; %stim as presented

    %convert to contrast (relative to mean) for filter convolution
    imMean = (luminanceData(ss).ImageMean  ./ luminanceData(ss).ImageMax);
    sStim = (sStim - imMean) / imMean;
    
    linearPrediction = conv(sStim,surround.LinearFilter);
    linearPrediction = linearPrediction(1:length(sStim));
    surroundGenSignal = cat(2,surroundGenSignal,linearPrediction);
end


figure(6); clf;
histogram2(centerGenSignal,surroundGenSignal,sqrt(numberOfBins_em),'DisplayStyle','tile',...
    'Normalization','probability','ShowEmptyBins','on');
% set(gca, 'ZScale', 'log')
% axis square; view(-5,70)
% xlim([-500 500]); ylim([-25 25])
xlabel('Center gen. signal'); ylabel('Surround gen. signal');
% zlabel('Prob');
%% look at errors of additivity for different models:
%SHARED NLINEARITY:
cc = linspace(min(centerGS),max(centerGS),sqrt(numberOfBins));
ss = linspace(min(surroundGS),max(surroundGS),sqrt(numberOfBins));
[CC,SS] = meshgrid(cc',ss');

resp_CenterOnly = CSModel_SharedNL(CC(:)',zeros(1,length(CC(:))),...
    fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
    fitRes_shared.gamma,fitRes_shared.epsilon);

resp_SurroundOnly = CSModel_SharedNL(zeros(1,length(CC(:))),SS(:)',...
    fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
    fitRes_shared.gamma,fitRes_shared.epsilon);

resp_CenterSurround = CSModel_SharedNL(CC(:)',SS(:)',...
    fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
    fitRes_shared.gamma,fitRes_shared.epsilon);

linearSummedResponse = resp_CenterOnly + resp_SurroundOnly;

err = resp_CenterSurround - linearSummedResponse;
errorMatrix = reshape(err,sqrt(numberOfBins),sqrt(numberOfBins));

figure(6); clf;
subplot(321)
plot(center.filterTimeVector,center.LinearFilter,'b')
hold on; plot([0 0.8],[0 0],'k:')
xlabel('Time (s)'); title('Center')
xlim([0 0.6]); ylim([min(center.LinearFilter) max(center.LinearFilter)])

subplot(322)
plot(surround.filterTimeVector,surround.LinearFilter,'r')
hold on; plot([0 0.8],[0 0],'k:')
xlabel('Time (s)'); title('Surround')
xlim([0 0.6]); ylim([min(surround.LinearFilter) max(surround.LinearFilter)])

subplot(3,2,3:6)
title('Shared model')
surfc(centerGS,surroundGS,errorMatrix);
% axis square; view(-16,22)
axis square; view(164,32)
xlabel('Center gen. signal'); ylabel('Surround gen. signal');
zlabel('Measured - Linear sum');

%JOINT NLINEARITY:
cc = linspace(min(centerGS),max(centerGS),sqrt(numberOfBins));
ss = linspace(min(surroundGS),max(surroundGS),sqrt(numberOfBins));
[CC,SS] = meshgrid(cc',ss');


resp_CenterOnly = JointNLin_mvcn(CC(:)',zeros(1,length(CC(:))),fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);

resp_SurroundOnly = JointNLin_mvcn(zeros(1,length(SS(:))),SS(:)',fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);

resp_CenterSurround = JointNLin_mvcn(CC(:)',SS(:)',fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);

linearSummedResponse = resp_CenterOnly + resp_SurroundOnly;

err = resp_CenterSurround - linearSummedResponse;
errorMatrix = reshape(err,sqrt(numberOfBins),sqrt(numberOfBins));

figure(7); clf;
subplot(321)
plot(center.filterTimeVector,center.LinearFilter,'b')
hold on; plot([0 0.8],[0 0],'k:')
xlabel('Time (s)'); title('Center')
xlim([0 0.6]); ylim([min(center.LinearFilter) max(center.LinearFilter)])

subplot(322)
title('Indep model')
plot(surround.filterTimeVector,surround.LinearFilter,'r')
hold on; plot([0 0.8],[0 0],'k:')
xlabel('Time (s)'); title('Surround')
xlim([0 0.6]); ylim([min(surround.LinearFilter) max(surround.LinearFilter)])

subplot(3,2,3:6)
surfc(centerGS,surroundGS,errorMatrix);
% axis square; view(-16,22)
axis square; view(164,32)
xlabel('Center gen. signal'); ylabel('Surround gen. signal');
zlabel('Measured - Linear sum');

%% all data - sandwich model:
load('horCRF.mat')
fitX = -1:0.01:1;
normFactor = CRFcumGauss(1,fitRes.alphaScale,fitRes.betaSens,fitRes.gammaXoff,fitRes.epsilonYoff);
fitY = CRFcumGauss(fitX,fitRes.alphaScale,fitRes.betaSens,fitRes.gammaXoff,fitRes.epsilonYoff) / normFactor;
% figure(5); clf; plot(fitX,fitY,'k-')

centerGS = center.generatorSignal(testDataInds);
surroundGS = surround.generatorSignal(testDataInds);
surroundGS = surroundGS ./ max(surroundGS);
measuredResponse = centerSurround.measuredResponse(testDataInds);
 
centerAloneResponse = center.measuredResponse(testDataInds);
surroundAloneResponse = surround.measuredResponse(testDataInds);
linearSummedResponse = centerAloneResponse + surroundAloneResponse;

aScale = -1;

upstreamSignal = centerGS - ...
    max(centerGS) * aScale .* CRFcumGauss(surroundGS,fitRes.alphaScale,fitRes.betaSens,fitRes.gammaXoff,fitRes.epsilonYoff) ./ normFactor;


sandwichResponse = sigmoidCRF(upstreamSignal,center.nonlinearity.fitParams.k,...
    center.nonlinearity.fitParams.c0,...
    center.nonlinearity.fitParams.amp,...
    center.nonlinearity.fitParams.yOff);

figure(6); clf;
subplot(311); hold on;
plot(measuredResponse,'k-')
plot(linearSummedResponse,'g')
ylim([min(measuredResponse) max(measuredResponse)])

subplot(312); hold on;
plot(measuredResponse,'k-')
plot(predictedResponse_shared,'b')
ylim([min(measuredResponse) max(measuredResponse)])

subplot(313); hold on;
plot(measuredResponse,'k-')
plot(sandwichResponse,'b')
ylim([min(measuredResponse) max(measuredResponse)])
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
PSTHsigma = 50; %msec
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

%% flag imageID node for example ON and OFF 

parentNode = gui.getSelectedEpochTreeNodes{1};

centerStims = {'Image','none','Image','Equiv','none','Equiv','Image','Equiv'};
surroundStims = {'none','Image','Image','none','Equiv','Equiv','Equiv','Image'};

populationNodes = {};
ct = 0;
for nn = 1:parentNode.descendentsDepthFirst.length
    if strcmp(parentNode.descendentsDepthFirst(nn).splitKey,...
            'protocolSettings(currentCenter)') && parentNode.descendentsDepthFirst(nn).custom.get('isSelected')
        ct = ct + 1;
        populationNodes(ct) = parentNode.descendentsDepthFirst(nn);
    end
end

clear responseMatrix
responseMatrix.ON = [];
responseMatrix.OFF = [];
for pp = 1:2
    if pp == 1
        fieldName = 'OFF'; 
    elseif pp == 2
        fieldName = 'ON'; 
    end
    for cc = 1:8
        currentNode = populationNodes{pp}.childBySplitValue(centerStims{cc}).childBySplitValue(surroundStims{cc});
        for ii = 1:currentNode.children.length
            newResp = getResponseAmplitudeStats(currentNode.children(ii).epochList,'extracellular');
            responseMatrix.(fieldName)(ii,cc) = newResp.integrated.mean;
        end
    end
end
%% Figures for COSYNE abstract:

%Linearity of surround:
figure(2); clf;
fig1=gca;
set(fig1,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Center + image surround (spikes)')
set(get(fig1,'YLabel'),'String','Center + linear surround (spikes)')

colors = pmkmp(3);

addLineToAxis(responseMatrix.ON(:,3),responseMatrix.ON(:,7),...
    'ONdata',fig1,colors(1,:),'none','o')
addLineToAxis(responseMatrix.OFF(:,3),responseMatrix.OFF(:,7),...
    'OFFdata',fig1,colors(2,:),'none','o')

addLineToAxis([0 30],[0 30],...
    'unity',fig1,'k',':','none')

makeAxisStruct(fig1,'LECS_parasols' ,'COSYNE2017Figs')


%%
clc;
load('NaturalImageFlashLibrary_101716.mat')
imageNames = fieldnames(imageData);
resourcesDir = '~/Documents/MATLAB/turner-package/resources';
patchSize = 100; %pixels
MicronsPerPixel = 6.6;
patchesPerImage = 100;
ImageID = 'imk03093';

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
figure(2); clf;
subplot(211); hold on
plot(sigmas.Center,responses.Center,'b-o');
plot(sigmas.Surround,responses.Surround,'r-o');
subplot(212)
surf(sigmas.Center,sigmas.Surround,responses.CenterSurround)
xlabel('Center'); ylabel('Surround')





