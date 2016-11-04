%% Make receptive field components:
% Center and surround are combined linearly then rectified
% 1) Linear center, linear surround. Sum C + S then rectify (DoG LN model)
% 2) Nonlinear center, linear surround. Sum C + S then rectify
% 3) Linear center, nonlinear surround. Sum C + S then rectify
% 4) Nonlinear center, nonlinear surorund. Sum C + S then rectify

clear all; close all; clc
scaleFactor = 10; %um per arcmin
FilterSize = round(800 ./ scaleFactor); % microns -> arcmins

% RF properties
%stdevs, microns to arcmin... 
CenterSubunitSigma = round(12 ./ scaleFactor);
SurroundSubunitSigma = round(50 ./ scaleFactor);
CenterSigma = round(50 ./ scaleFactor);
SurroundSigma = round(150 ./ scaleFactor);

CenterWeight = 0.6;

% Center subunit locations - square grid
CenterSubunitFilter = zeros(FilterSize, FilterSize);
SubunitLocations = find(rem([1:FilterSize], 2*CenterSubunitSigma) == 0);
for x = 1:length(SubunitLocations)
    CenterSubunitFilter(SubunitLocations(x), SubunitLocations) = 1;
end
CenterSubunitIndices = find(CenterSubunitFilter > 0);

% Surround subunit locations - square grid
SurroundSubunitFilter = zeros(FilterSize, FilterSize);
SubunitLocations = find(rem([1:FilterSize], 2*SurroundSubunitSigma) == 0);
for x = 1:length(SubunitLocations)
    SurroundSubunitFilter(SubunitLocations(x), SubunitLocations) = 1;
end
SurroundSubunitIndices = find(SurroundSubunitFilter > 0);


% Filters for: 
% 1) center subunit
% 2) surround subunit
% 3) linear center
% 4) linear surround
for x = 1:FilterSize
    for y = 1:FilterSize
        CenterSubunit(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (CenterSubunitSigma^2)));
        SurroundSubunit(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (SurroundSubunitSigma^2)));
        Center(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (CenterSigma^2)));
        Surround(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (SurroundSigma^2)));
    end
end
centerWeights = Center(CenterSubunitIndices);
surroundWeights = Surround(SurroundSubunitIndices);

% normalize each component
centerWeights = centerWeights / sum(centerWeights);
surroundWeights = surroundWeights / sum(surroundWeights);

CenterSubunit = CenterSubunit / sum(CenterSubunit(:));
SurroundSubunit = SurroundSubunit / sum(SurroundSubunit(:));
Center = Center / sum(Center(:));
Surround = Surround / sum(Surround(:));
% 
% figure(1); clf;
% subplot(221); imagesc(CenterSubunit); title('CenterSubunit')
% subplot(222); imagesc(SurroundSubunit);  title('SurroundSubunit')
% subplot(223); imagesc(Center); title('Linear center')
% subplot(224); imagesc(Surround); title('Linear surround')

%% Expanding spots:

spotDiameter = round((50:50:800) ./ scaleFactor);
spotIntensity = 1;
[rr, cc] = meshgrid(1:FilterSize,1:FilterSize);

response.LinearCenter = [];
response.NonlinearCenter = [];
response.LinearSurround = [];
response.NonlinearSurround = [];
for ss = 1:length(spotDiameter)
    currentRadius = spotDiameter(ss)/2;
    stimulusImage = double(sqrt((rr-round(FilterSize/2)).^2+(cc-round(FilterSize/2)).^2)<=currentRadius);
    % CENTER
    % Convolve with center subunits
    centerSubunitFilteredImage = conv2(stimulusImage, CenterSubunit, 'same');
    centerSubunitActivations = centerSubunitFilteredImage(CenterSubunitIndices);
    % Linear center:
    response.LinearCenter(ss) = sum(centerSubunitActivations .* centerWeights);
    % Nonlinear center:
    centerSubunitOutputs = centerSubunitActivations;
    centerSubunitOutputs(centerSubunitOutputs < 0) = 0;
    response.NonlinearCenter(ss) = sum(centerSubunitOutputs .* centerWeights);
    
    % SURROUND
    % Convolve with surround subunits
    surroundSubunitFilteredImage = conv2(stimulusImage, SurroundSubunit, 'same');
    surroundSubunitActivations = surroundSubunitFilteredImage(SurroundSubunitIndices);
    % Linear surround:
    response.LinearSurround(ss) = sum(surroundSubunitActivations .* surroundWeights);
    % Nonlinear surround:
    surroundSubunitOutputs = surroundSubunitActivations;
    surroundSubunitOutputs(surroundSubunitOutputs < 0) = 0;
    response.NonlinearSurround(ss) = sum(surroundSubunitOutputs .* surroundWeights);
end

spotDiameter = spotDiameter .* scaleFactor;
figure(2); clf; hold on;
LL = CenterWeight .* response.LinearCenter - (1 - CenterWeight) .* response.LinearSurround;
LN = CenterWeight .* response.LinearCenter - (1 - CenterWeight) .* response.NonlinearSurround;
NL = CenterWeight .* response.NonlinearCenter - (1 - CenterWeight) .* response.LinearSurround;
NN = CenterWeight .* response.NonlinearCenter - (1 - CenterWeight) .* response.NonlinearSurround;

plot(spotDiameter, LL, 'ko');
plot(spotDiameter, LN, 'ko');
plot(spotDiameter, NL, 'ko');
plot(spotDiameter, NN, 'ko');

%% Contrast reversing gratings:
barWidth = round((20:20:400) ./ scaleFactor);
gratingContrast = 0.9;
[rr, cc] = meshgrid(1:FilterSize,1:FilterSize);

response.LinearCenter = [];
response.NonlinearCenter = [];
response.LinearSurround = [];
response.NonlinearSurround = [];
for bb = 1:length(barWidth)
    tempStim_right = sin(2*(1/(2*barWidth(bb)))*pi.*(0:FilterSize/2-1));
    tempStim_left = sin(2*(1/(2*barWidth(bb)))*pi.*(1:FilterSize/2));
    
    tempStim = [-fliplr(tempStim_left) tempStim_right];
    
    tempStim(tempStim>0) = gratingContrast; tempStim(tempStim<=0) = -gratingContrast; %contrast
    stimulusImage = repmat(tempStim,FilterSize,1);
    
    % CENTER
    % Convolve with center subunits
    centerSubunitFilteredImage = conv2(stimulusImage, CenterSubunit, 'same');
    centerSubunitActivations = centerSubunitFilteredImage(CenterSubunitIndices);
    % Linear center:
    response.LinearCenter(bb) = sum(centerSubunitActivations .* centerWeights);
    % Nonlinear center:
    centerSubunitOutputs = centerSubunitActivations;
    centerSubunitOutputs(centerSubunitOutputs < 0) = 0;
    response.NonlinearCenter(bb) = sum(centerSubunitOutputs .* centerWeights);
    
    % SURROUND
    % Convolve with surround subunits
    surroundSubunitFilteredImage = conv2(stimulusImage, SurroundSubunit, 'same');
    surroundSubunitActivations = surroundSubunitFilteredImage(SurroundSubunitIndices);
    % Linear surround:
    response.LinearSurround(bb) = sum(surroundSubunitActivations .* surroundWeights);
    % Nonlinear surround:
    surroundSubunitOutputs = surroundSubunitActivations;
    surroundSubunitOutputs(surroundSubunitOutputs < 0) = 0;
    response.NonlinearSurround(bb) = sum(surroundSubunitOutputs .* surroundWeights);
end

barWidth = barWidth .* scaleFactor;
figure(3); clf; 
subplot(211); plot(barWidth, response.LinearCenter, 'ko');
hold on; plot(barWidth, response.NonlinearCenter, 'rx');
title('Center')
xlabel('Bar width (um'); ylabel('Response'); legend('Linear','Nonlinear')

subplot(212); plot(barWidth, response.LinearSurround, 'ko');
hold on; plot(barWidth, response.NonlinearSurround, 'rx');
title('Surround')
xlabel('Bar width (um'); ylabel('Response'); legend('Linear','Nonlinear')

%% Natural image patches:
ImageIndex = 7;
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
figure(4); 
clf;  
imagesc((my_image.^0.3)');colormap gray;axis image; axis off; hold on;
[ImageX, ImageY] = size(my_image);

% scale image to [0 1] -- contrast, relative to mean over entire image...
my_image_nomean = (my_image - mean(my_image(:))) ./ mean(my_image(:));

% measure RF component responses to patches
patchLocation = [];
response.LinearCenter = [];
response.NonlinearCenter = [];
response.LinearSurround = [];
response.NonlinearSurround = [];
for patch = 1:NumFixations
    % choose location
    x = round(FilterSize/2 + (ImageX - FilterSize)*rand);
    y = round(FilterSize/2 + (ImageY - FilterSize)*rand);
    patchLocation(patch) = x;
    patchLocation(patch) = y;
    stimulusImage = my_image_nomean(x-round(FilterSize/2)+1:x+round(FilterSize/2),y-round(FilterSize/2)+1:y+round(FilterSize/2));
    
    % CENTER
    % Convolve with center subunits
    centerSubunitFilteredImage = conv2(stimulusImage, CenterSubunit, 'same');
    centerSubunitActivations = centerSubunitFilteredImage(CenterSubunitIndices);
    % Linear center:
    response.LinearCenter(patch) = sum(centerSubunitActivations .* centerWeights);
    % Nonlinear center:
    centerSubunitOutputs = centerSubunitActivations;
    centerSubunitOutputs(centerSubunitOutputs < 0) = 0;
    response.NonlinearCenter(patch) = sum(centerSubunitOutputs .* centerWeights);
    
    % SURROUND
    % Convolve with surround subunits
    surroundSubunitFilteredImage = -conv2(stimulusImage, SurroundSubunit, 'same');
    surroundSubunitActivations = surroundSubunitFilteredImage(SurroundSubunitIndices);
    % Linear surround:
    response.LinearSurround(patch) = sum(surroundSubunitActivations .* surroundWeights);
    % Nonlinear surround:
    surroundSubunitOutputs = surroundSubunitActivations;
    surroundSubunitOutputs(surroundSubunitOutputs < 0) = 0;
    response.NonlinearSurround(patch) = sum(surroundSubunitOutputs .* surroundWeights);
    
    
    if (rem(patch, 500) == 0)
        fprintf(1, '%d ', patch);
    end
end
response.LinearCenter = max(response.LinearCenter,0);
response.LinearSurround = max(response.LinearSurround,0);
%% plot image patch responses

% figure(5); clf; 
% subplot(121); plot(response.LinearCenter, response.NonlinearCenter, 'ko');
% xlabel('Linear Center'); ylabel('Nonlinear Center')
% subplot(122); plot(response.LinearSurround, response.NonlinearSurround, 'ko');
% xlabel('Linear Surround'); ylabel('Nonlinear Surround')

LL = CenterWeight .* response.LinearCenter + (1 - CenterWeight) .* response.LinearSurround;
LN = CenterWeight .* response.LinearCenter + (1 - CenterWeight) .* response.NonlinearSurround;
NL = CenterWeight .* response.NonlinearCenter + (1 - CenterWeight) .* response.LinearSurround;
NN = CenterWeight .* response.NonlinearCenter + (1 - CenterWeight) .* response.NonlinearSurround;

% LL = max(LL,0); LN = max(LN,0); NL = max(NL,0); NN = max(NN,0);

figure(6); clf;
subplot(141); plot(max(response.NonlinearCenter,0),...
    max(response.LinearCenter,0), 'ko');
xlabel('Nonlinear center'); ylabel('Linear center')

subplot(142); plot(NL, LL, 'ko');
xlabel('Nonlinear Center + linear surround'); ylabel('Linear Center + linear surround')
subplot(143); plot(NN, LL, 'ko');
xlabel('Nonlinear Center + nonlinear surround'); ylabel('Linear Center + linear surround')

subplot(144); plot(LN, LL, 'ko');
xlabel('Linear Center + nonlinear surround'); ylabel('Linear Center + linear surround')



figure(7); clf;
subplot(121); plot(max(response.NonlinearSurround,0),...
    max(response.LinearSurround,0), 'ko');
xlabel('Nonlinear surround'); ylabel('Linear surround')

subplot(122); plot(LN, LL, 'ko');
xlabel('Linear Center + nonlinear surround'); ylabel('Linear Center + linear surround')

% response histograms:
yUP = 140;
figure(9); clf;
subplot(411)
hold on
tempCenter = max(response.LinearCenter,0);
tempCenter = tempCenter ./ max(tempCenter);
[n edges] = histcounts(tempCenter,40);
plot(edges(2:end),n,'k');
tempCS = LL;
tempCS = tempCS ./ max(tempCS);
[n edges] = histcounts(tempCS,40);
plot(edges(2:end),n,'r');
ylim([0 yUP])
legend('Linear center','Linear center + linear surround')

subplot(412)
hold on
tempCenter = max(response.LinearCenter,0);
tempCenter = tempCenter ./ max(tempCenter);
[n edges] = histcounts(tempCenter,40);
plot(edges(2:end),n,'k');
tempCS = LN;
tempCS = tempCS ./ max(tempCS);
[n edges] = histcounts(tempCS,40);
plot(edges(2:end),n,'r');
ylim([0 yUP])
legend('Linear center','Linear center + nonlinear surround')

subplot(413)
hold on
tempCenter = max(response.NonlinearCenter,0);
tempCenter = tempCenter ./ max(tempCenter);
[n edges] = histcounts(tempCenter,40);
plot(edges(2:end),n,'k');
tempCS = NL;
tempCS = tempCS ./ max(tempCS);
[n edges] = histcounts(tempCS,40);
plot(edges(2:end),n,'r');
ylim([0 yUP])
legend('Nonlinear center','Nonlinear center + linear surround')

subplot(414)
hold on
tempCenter = max(response.NonlinearCenter,0);
tempCenter = tempCenter ./ max(tempCenter);
[n edges] = histcounts(tempCenter,40);
plot(edges(2:end),n,'k');
tempCS = NN;
tempCS = tempCS ./ max(tempCS);
[n edges] = histcounts(tempCS,40);
plot(edges(2:end),n,'r');
ylim([0 yUP])
legend('Nonlinear center','Nonlinear center + Nonlinear surround')

%% corrs

intensityProxy = response.LinearCenter;

figure(10); clf;
subplot(511);
plot(intensityProxy, max(response.NonlinearCenter,0),'ko')
xlabel('Center Intensity'); ylabel('Resp')
title('Nonlinear center')

subplot(512);
plot(intensityProxy, LL,'ko')
xlabel('Center Intensity'); ylabel('Resp')
title('Linear center & surround')

subplot(513);
plot(intensityProxy, NL,'ko')
xlabel('Center Intensity'); ylabel('Resp')
title('Nonlinear center & linear surround')

subplot(514);
plot(intensityProxy, LN,'ko')
xlabel('Center Intensity'); ylabel('Resp')
title('Linear center & nonlinear surround')

subplot(515);
plot(intensityProxy, NN,'ko')
xlabel('Center Intensity'); ylabel('Resp')
title('Nonlinear center & surround')

%% mean + grating to find subunit polarity

meanIntensity = -0.9:0.1:0.9;
barWidth = round(200./ scaleFactor);
gratingContrast = 0.9;
[rr, cc] = meshgrid(1:FilterSize,1:FilterSize);

response.OFFsubunits = [];
response.ONsubunits = [];
response.LinearSurround = [];
for bb = 1:length(meanIntensity)
    tempStim_right = sin(2*(1/(2*barWidth))*pi.*(0:FilterSize/2-1));
    tempStim_left = sin(2*(1/(2*barWidth))*pi.*(1:FilterSize/2));
    
    tempStim = [-fliplr(tempStim_left) tempStim_right];
    
    tempStim(tempStim>0) = meanIntensity(bb) + gratingContrast;
    tempStim(tempStim<=0) = meanIntensity(bb) - gratingContrast; %contrast
    stimulusImage = repmat(tempStim,FilterSize,1);

    
    % SURROUND
    % Convolve with surround subunits
    surroundSubunitFilteredImage = conv2(stimulusImage, SurroundSubunit, 'same');
    surroundSubunitActivations = surroundSubunitFilteredImage(SurroundSubunitIndices);
    
    % Linear (ON) surround:
    response.LinearSurround(bb) = sum(surroundSubunitActivations .* surroundWeights);
    
    % ON subunits:
    surroundSubunitOutputs = surroundSubunitActivations;
    surroundSubunitOutputs(surroundSubunitOutputs < 0) = 0;
    response.ONsubunits(bb) = sum(surroundSubunitOutputs .* surroundWeights);
    % OFF subunits:
    surroundSubunitOutputs = -surroundSubunitActivations;
    surroundSubunitOutputs(surroundSubunitOutputs < 0) = 0;
    response.OFFsubunits(bb) = sum(surroundSubunitOutputs .* surroundWeights);
end

figure(1); clf; hold on;
plot(meanIntensity, max(response.LinearSurround,0),'ko')
plot(meanIntensity, max(response.LinearSurround + response.ONsubunits,0),'bx')
plot(meanIntensity, max(response.LinearSurround + response.OFFsubunits,0),'rx')
legend('Offset alone','Offset + grating (ON)', 'Offset + grating (OFF)')
xlabel('Offset Intensity'); ylabel('Response');


