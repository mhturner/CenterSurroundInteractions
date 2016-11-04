%% Make receptive field components:
% two models:
% 1) Indep nonlinearity, f(subunit center + linear surround)
% 2) Shared nonlinearity, f(DoG subunit)
% 3) LN model (DoG)
% 4) Linear center only
% 5) Subunit center only
clear all; close all; clc
scaleFactor = 10; %um per arcmin
FilterSize = round(800 ./ scaleFactor); % microns -> arcmins

% RF properties
%stdevs, microns to arcmin... 
SubunitCenterSigma = round(6 ./ scaleFactor);
SubunitSurroundSigma = round(150 ./ scaleFactor);
LinearCenterSigma = round(40 ./ scaleFactor);
LinearSurroundSigma = round(150 ./ scaleFactor);

SubunitCenterWeight = 0.6;
LinearCenterWeight = 0.73;

% bipolar subunit locations - square grid
TempFilter = zeros(FilterSize, FilterSize);
SubunitLocations = find(rem([1:FilterSize], 2*SubunitCenterSigma) == 0);
for x = 1:length(SubunitLocations)
    TempFilter(SubunitLocations(x), SubunitLocations) = 1;
end
SubunitIndices = find(TempFilter > 0);

% Filters for: 
% 1) subunit (bipolar) center
% 2) subunit (bipolar) surround
% 3) linear center
% 4) linear surround
for x = 1:FilterSize
    for y = 1:FilterSize
        SubunitCenter(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (SubunitCenterSigma^2)));
        SubunitSurround(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (SubunitSurroundSigma^2)));
        LinearCenter(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (LinearCenterSigma^2)));
        LinearSurround(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (LinearSurroundSigma^2)));
    end
end
subunitWeights = LinearCenter(SubunitIndices);

% normalize each component
subunitWeights = subunitWeights / sum(subunitWeights);

SubunitCenter = SubunitCenterWeight .* SubunitCenter / sum(SubunitCenter(:));
SubunitSurround = -(1 - SubunitCenterWeight) * SubunitSurround / sum(SubunitSurround(:));
LinearCenter = LinearCenterWeight * LinearCenter / sum(LinearCenter(:));
LinearSurround = -(1 - LinearCenterWeight) * LinearSurround / sum(LinearSurround(:));

Subunit_DoG = SubunitCenterWeight * SubunitCenter + (1 - SubunitCenterWeight) * SubunitSurround;

figure(1); clf;
subplot(231); imagesc(SubunitCenter); title('Subunit center')
subplot(232); imagesc(SubunitSurround);  title('Subunit surround')
subplot(233); imagesc(Subunit_DoG);  title('Subunit DoG')

subplot(234); imagesc(LinearCenter); title('Linear center')
subplot(235); imagesc(LinearSurround); title('Linear surround')
subplot(236); imagesc(LinearCenter + LinearSurround); title('Linear DoG')

%% Expanding spots:
spotDiameter = round((50:50:800) ./ scaleFactor);
spotIntensity = 1;
[rr, cc] = meshgrid(1:FilterSize,1:FilterSize);
response.IndepNonlinearity = [];
response.SharedNonlinearity = [];
for ss = 1:length(spotDiameter)
    currentRadius = spotDiameter(ss)/2;
    stimulusImage = double(sqrt((rr-round(FilterSize/2)).^2+(cc-round(FilterSize/2)).^2)<=currentRadius);
    % Get activation of linear surround
    LinearSurroundResponse = sum(sum(stimulusImage .* LinearSurround));
    
    % get activation of each subunit center:
    subunitFilteredImage = conv2(stimulusImage, SubunitCenter, 'same');
    subunitCenterActivations = subunitFilteredImage(SubunitIndices);
    
    % get activation of each subunit DoG:
    subunitFilteredImage = conv2(stimulusImage, Subunit_DoG, 'same');
    subunitDoGActivations = subunitFilteredImage(SubunitIndices);
    
    %model 1: nonlinear center, linear surround. Rectify center subunits
    %and sum, then add surround to center.
    subunitOutputs = subunitCenterActivations;
    subunitOutputs(subunitOutputs<0) = 0;
    temp_Center = sum(subunitOutputs .* subunitWeights);
    temp_Surround = LinearSurroundResponse;
    response.IndepNonlinearity(ss) = max(temp_Center + temp_Surround,0); %rectify output

    %model 2: nonlinear subunits with a DoG RF. Rectify then sum
    subunitOutputs = subunitDoGActivations;
    subunitOutputs(subunitOutputs < 0) = 0;
    response.SharedNonlinearity(ss) = sum(subunitOutputs.* subunitWeights);
end

spotDiameter = spotDiameter .* scaleFactor;
figure(2); clf; hold on;
plot(spotDiameter, response.IndepNonlinearity ./ max(response.IndepNonlinearity), 'ko')
plot(spotDiameter, response.SharedNonlinearity ./ max(response.SharedNonlinearity), 'rx')
xlabel('Spot diameter (um)'); ylabel('Resp (normalized)');
legend('Ind. NL','Shared NL (DoG subunits)')

%% Contrast reversing gratings:
barWidth = round((20:20:400) ./ scaleFactor);
gratingContrast = 0.9;
[rr, cc] = meshgrid(1:FilterSize,1:FilterSize);

response.IndepNonlinearity = [];
response.SharedNonlinearity = [];
for bb = 1:length(barWidth)
    tempStim_right = sin(2*(1/(2*barWidth(bb)))*pi.*(0:FilterSize/2-1));
    tempStim_left = sin(2*(1/(2*barWidth(bb)))*pi.*(1:FilterSize/2));
    
    tempStim = [-fliplr(tempStim_left) tempStim_right];
    
    tempStim(tempStim>0) = gratingContrast; tempStim(tempStim<=0) = -gratingContrast; %contrast
    stimulusImage = repmat(tempStim,FilterSize,1);
    
    % Get activation of linear surround
    LinearSurroundResponse = sum(sum(stimulusImage .* LinearSurround));
    
    % get activation of each subunit center:
    subunitFilteredImage = conv2(stimulusImage, SubunitCenter, 'same');
    subunitCenterActivations = subunitFilteredImage(SubunitIndices);
    
    % get activation of each subunit DoG:
    subunitFilteredImage = conv2(stimulusImage, Subunit_DoG, 'same');
    subunitDoGActivations = subunitFilteredImage(SubunitIndices);
    
    %model 1: nonlinear center, linear surround. Rectify center subunits
    %and sum, then add surround to center.
    subunitOutputs = subunitCenterActivations;
    subunitOutputs(subunitOutputs<0) = 0;
    temp_Center = sum(subunitOutputs .* subunitWeights);
    temp_Surround = LinearSurroundResponse;
    response.IndepNonlinearity(bb) = max(temp_Center + temp_Surround,0); %rectify output

    %model 2: nonlinear subunits with a DoG RF. Rectify then sum
    subunitOutputs = subunitDoGActivations;
    subunitOutputs(subunitOutputs < 0) = 0;
    response.SharedNonlinearity(bb) = sum(subunitOutputs.* subunitWeights);
end

barWidth = barWidth .* scaleFactor;
figure(3); clf;
plot(barWidth, response.IndepNonlinearity ./ max(response.IndepNonlinearity), 'ko');
hold on; plot(barWidth, response.SharedNonlinearity ./ max(response.SharedNonlinearity), 'rx');
title('Center')
xlabel('Bar width (um'); ylabel('Response');
legend('Ind. NL','Shared NL (DoG subunits)')


%% Natural image patches:
ImageIndex = 1; % 3
NumFixations = 1000;           % how many patches to sample
randSeed = 1;
rng(randSeed);

% Pull van Hateren images...
IMAGES_DIR            = '~/Documents/MATLAB/Analysis/NatImages/Stimuli/VHsubsample_20160105/';
temp_names                  = GetFilenames(IMAGES_DIR,'.iml');
for file_num = 1:size(temp_names,1)
    temp                    = temp_names(file_num,:);
    temp                    = deblank(temp);
    img_filenames_list{file_num}  = temp;
end
img_filenames_list = sort(img_filenames_list);
clc;

% Load  and plot the image to analyze
f1=fopen([IMAGES_DIR, img_filenames_list{ImageIndex}],'rb','ieee-be');
w=1536;h=1024;
my_image=fread(f1,[w,h],'uint16');
figure(6); 
clf;  
imagesc((my_image.^0.3)');colormap gray;axis image; axis off; hold on;
[ImageX, ImageY] = size(my_image);

% scale image to [0 1] -- contrast, relative to mean over entire image...
my_image_nomean = (my_image - mean(my_image(:))) ./ mean(my_image(:));

% measure RF component responses to patches
patchLocation = [];
response.IndepNonlinearity = [];
response.SharedNonlinearity = [];
response.DoG = [];
response.LinearCenter = [];
response.SubunitCenter = [];
response.LNCenterPlusSurround = [];
for patch = 1:NumFixations
    % choose location
    x = round(FilterSize/2 + (ImageX - FilterSize)*rand);
    y = round(FilterSize/2 + (ImageY - FilterSize)*rand);
    patchLocation(1,patch) = x;
    patchLocation(2,patch) = y;
    stimulusImage = my_image_nomean(x-round(FilterSize/2)+1:x+round(FilterSize/2),y-round(FilterSize/2)+1:y+round(FilterSize/2));
    
    % Get activation of linear surround
    LinearSurroundResponse = sum(sum(stimulusImage .* LinearSurround));
    
    % get activation of each subunit center:
    subunitFilteredImage = conv2(stimulusImage, SubunitCenter, 'same');
    subunitCenterActivations = subunitFilteredImage(SubunitIndices);
    
    % get activation of each subunit DoG:
    subunitFilteredImage = conv2(stimulusImage, Subunit_DoG, 'same');
    subunitDoGActivations = subunitFilteredImage(SubunitIndices);
    
    %model 1: nonlinear center, linear surround. Rectify center subunits
    %and sum, then add surround to center.
    subunitOutputs = subunitCenterActivations;
    response.LinearCenter(patch) = sum(subunitOutputs .* subunitWeights);
    response.LNCenterPlusSurround(patch) = ...
        max(max(response.LinearCenter(patch),0) + LinearSurroundResponse,0);
    %rectify for nonlinear center: 
    subunitOutputs(subunitOutputs<0) = 0; 
    temp_Center = sum(subunitOutputs .* subunitWeights);
    response.SubunitCenter(patch) = temp_Center;
    temp_Surround = LinearSurroundResponse;
    response.IndepNonlinearity(patch) = max(temp_Center + temp_Surround,0); %rectify output
    
    

    %model 2: nonlinear subunits with a DoG RF. Rectify then sum
    subunitOutputs = subunitDoGActivations;
    subunitOutputs(subunitOutputs < 0) = 0;
    response.SharedNonlinearity(patch) = sum(subunitOutputs.* subunitWeights);
    
    %model 3: LN model, DoG
    temp_Center = sum(subunitCenterActivations .* subunitWeights);
    response.DoG(patch) = max(temp_Center - LinearSurroundResponse,0);
    
    if (rem(patch, 500) == 0)
        fprintf(1, '%d ', patch);
    end
end


%%
figure(4); clf;
subplot(221);
plot(response.SubunitCenter, max(response.LinearCenter,0),'ko')
xlabel('Subunit center'); ylabel('Linear center')

subplot(222);
plot(response.IndepNonlinearity, response.DoG,'ko')
xlabel('Indep nonlinearity'); ylabel('Linear DoG')

subplot(223);
plot(response.SharedNonlinearity, response.DoG,'ko')
xlabel('Shared nonlinearity'); ylabel('Linear DoG')

subplot(224);
plot(response.IndepNonlinearity, response.SharedNonlinearity,'ko')
xlabel('Indep nonlinearity'); ylabel('Shared nonlinearity');



intensityProxy = response.LinearCenter;
figure(5); clf

subplot(421);
plot(intensityProxy, response.SubunitCenter,'ko')
title('SubunitCenter')
subplot(422)
[n, edges] = histcounts(response.SubunitCenter,40);
plot(edges(2:end),n,'k');
ylim([0 250])

subplot(423);
plot(intensityProxy, response.DoG,'ko')
title('DoG')
subplot(424)
[n, edges] = histcounts(response.DoG,40);
plot(edges(2:end),n,'k');
ylim([0 250])


subplot(425);
plot(intensityProxy, response.IndepNonlinearity,'ko')
title('IndepNonlinearity')
subplot(426)
[n, edges] = histcounts(response.IndepNonlinearity,40);
plot(edges(2:end),n,'k');
ylim([0 250])

subplot(427);
plot(intensityProxy, response.SharedNonlinearity,'ko')
xlabel('Intensity in center'); ylabel('Response')
title('SharedNonlinearity')
subplot(428)
[n, edges] = histcounts(response.SharedNonlinearity,40);
plot(edges(2:end),n,'k');
xlabel('Response'); ylabel('Count')
ylim([0 250])


%% diff heat maps
mapVals = response.IndepNonlinearity - response.SharedNonlinearity;

figure(6); 
clf;  
imagesc((my_image.^0.3)');colormap gray;axis image; axis off; hold on;
freezeColors
hold on;
scaleFactor = max(mapVals);
cColors =mapVals./scaleFactor;
scatter(patchLocation(1,:),patchLocation(2,:),15,cColors,'filled')
set(gca,'clim',[min(mapVals) max(cColors)])
colormap(jet); colorbar

