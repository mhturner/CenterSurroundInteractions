%% Compute new RF surround model results
clear all; close all; clc;

outputDir = '~/Dropbox/RiekeLab/Analysis/MATLAB/RFSurround/resources/RFsurroundNLIresults/';
resources_dir = '~/Dropbox/RiekeLab/Analysis/MATLAB/turner-package/resources/';
IMAGES_DIR_VH = '~/Dropbox/RiekeLab/Analysis/MATLAB/MHT-analysis/resources/vanhateren_iml/';
load([resources_dir, 'dovesFEMstims_20160422.mat'])

targetImageIndex = 1;
noPatches = 50;

micronsPerPixel = 6.6;

%RF components:
%"Standard" = 0.72
% Weak = 0.10
% Strong = 1.1
SubunitSurroundWeight = 0.72;     % relative to center integral

filterSize_um = 700;                % size of patch (um)
subunitSigma_um = 10;
subunitSurroundSigma_um = 150;
centerSigma_um = 40;

% % % % % % % % LOAD NATURAL IMAGE % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for ss = 1:length(FEMdata)
    if FEMdata(ss).ImageIndex == targetImageIndex
        stimInd = ss;
        break;
    end
end
%image:
fileName = [FEMdata(stimInd).ImageName(1:8), '_',num2str(FEMdata(stimInd).SubjectIndex),'.mat'];
f1=fopen([IMAGES_DIR_VH, FEMdata(stimInd).ImageName],'rb','ieee-be');
w=1536;h=1024;
my_image=fread(f1,[w,h],'uint16');
my_image = my_image';
my_image = my_image ./ max(my_image(:)); %normalize to [0 1]

my_image_nomean = (my_image - mean(my_image(:))) ./ mean(my_image(:));

%convert to pixels:
filterSize = round(filterSize_um / micronsPerPixel);                % size of patch (um)
subunitSigma = round(subunitSigma_um / micronsPerPixel);
subunitSurroundSigma = round(subunitSurroundSigma_um / micronsPerPixel);
centerSigma = round(centerSigma_um / micronsPerPixel);

% % % % % % % % SPATIAL RF % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
TempFilter = zeros(filterSize, filterSize);
SubunitLocations = find(rem([1:filterSize], 2*subunitSigma) == 0);
for x = 1:length(SubunitLocations)
    TempFilter(SubunitLocations(x), SubunitLocations) = 1;
end
SubunitIndices = find(TempFilter > 0);

% Filters for subunit with surround, center, surround
for x = 1:filterSize
    for y = 1:filterSize
        SubunitFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (subunitSigma^2)));
        SubunitSurroundFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (subunitSurroundSigma^2)));
        RFCenter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (centerSigma^2)));
    end
end
% normalize components
SubunitFilter = SubunitFilter / sum(SubunitFilter(:));
SubunitSurroundFilter = SubunitSurroundFilter / sum(SubunitSurroundFilter(:));
SubunitWithSurroundFilter = SubunitFilter - SubunitSurroundWeight .* SubunitSurroundFilter; %c/s subunit filter
RFCenter = RFCenter / sum(RFCenter(:));
%get weighting of each subunit output by center
subunitWeightings = RFCenter(SubunitIndices);
subunitWeightings = subunitWeightings ./ sum(subunitWeightings);

% % % % % % % % TEMPORAL RF % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
load('OffParasolExcitatoryFilters.mat')
tempNorm = filters.center ./ repmat(max(abs(filters.center),[],2),1,size(filters.center,2));
centerFilter = mean(tempNorm); centerFilter = -centerFilter ./ max(abs(centerFilter));
tempNorm = filters.surround ./ repmat(max(abs(filters.surround),[],2),1,size(filters.surround,2));
surroundFilter = mean(tempNorm); surroundFilter = surroundFilter ./ max(abs(surroundFilter));

figure(1); clf; hold on;
plot(centerFilter,'b');
plot(surroundFilter,'r');

% % % % % % % % RESPONSE TO SPOTS % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

stimulus.spotSize = 0:10:700; %diameter of spot

response.spots = [];

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
colors = pmkmp(length(stimulus.spotSize));
figure(2); clf; hold on;
for ss = 1:length(stimulus.spotSize) %get responses to each spot
    currentRadius = (stimulus.spotSize(ss)/2)/ micronsPerPixel; %convert to pixel
    spotBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=currentRadius);
    
    %activation in space:
    convolved_SubunitCenter = conv2(spotBinary, SubunitFilter, 'same');
    subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
    convolved_SubunitSurround = conv2(spotBinary, SubunitSurroundFilter, 'same');
    subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
    %temporal activation of each subunit center & surround
    %rows = subunit, cols = time
    subunitVoltage = subunitCenterActivations * centerFilter -...
        subunitSuroundActivations * surroundFilter;

    %hit each subunit with rectifying synapse
    rectInds = find(subunitVoltage < 0);
    subunitVoltage(rectInds) = 0;
    subunitOutputs = subunitVoltage;
    response.spots(ss,:) = sum(subunitOutputs .* subunitWeightings); %sum over subunits
    plot(response.spots(ss,:),'Color',colors(ss,:));
end

intResponse = trapz(response.spots(:,1:1000),2);
figure(3); clf; plot(stimulus.spotSize,intResponse,'ko')


[~, ind] = max(intResponse);
centerSize_um = stimulus.spotSize(ind);
centerSize = round(centerSize_um / micronsPerPixel);

%%
% % % responses to image patches and discs, +/- surrounds

tic;
contrastPolarity = -1;

rng(1); %set random seed

%center only:
response.Image = [];
response.Disc = [];
%matched surround:
response.ImageMatchedSurround = [];
response.DiscMatchedSurround = [];
%mixed surround:
response.ImageMixedSurround = [];
response.DiscMixedSurround = [];

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
centerBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=(centerSize/2));
surroundBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)>(centerSize/2));
for pp = 1:noPatches
    %main image patch location:
    x = round(filterSize/2 + (w - filterSize)*rand);
    y = round(filterSize/2 + (h - filterSize)*rand);
    newImagePatch = contrastPolarity .* my_image_nomean(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    stimulus.CenterLocation(pp,:) = [x, y];

    %mixed image surround location:
    x_mix = round(filterSize/2 + (w - filterSize)*rand);
    y_mix = round(filterSize/2 + (h - filterSize)*rand);
    mixedImagePatch = contrastPolarity .* my_image_nomean(y_mix - filterSize/2 + 1 : y_mix + filterSize/2, x_mix - filterSize/2 + 1 : x_mix + filterSize/2);
    stimulus.MixedSurroundLocation(pp,:) = [x_mix, y_mix];
    

    %Stim responses:
    %1) Image:
        ImageStim = newImagePatch;
        ImageStim(~centerBinary) = 0; %add aperture
        %activation in space:
        convolved_SubunitCenter = conv2(ImageStim, SubunitFilter, 'same');
        subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
        convolved_SubunitSurround = conv2(ImageStim, SubunitSurroundFilter, 'same');
        subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
        %temporal activation of each subunit center & surround
        %rows = subunit, cols = time
        subunitVoltage = subunitCenterActivations * centerFilter -...
            subunitSuroundActivations * surroundFilter;
        %hit each subunit with rectifying synapse
        rectInds = find(subunitVoltage < 0);
        subunitVoltage(rectInds) = 0;
        subunitOutputs = subunitVoltage;
        response.Image(pp,:) = sum(subunitOutputs .* subunitWeightings);

        

    %2) Disc:
        %compute linear equivalent contrast...
        % Get the model RF:
        RF = fspecial('gaussian',filterSize,centerSigma);
        weightingFxn = centerBinary .* RF; %set to zero mean gray pixels
        weightingFxn = weightingFxn ./ sum(weightingFxn(:)); %sum to one
        equivalentContrast = sum(sum(weightingFxn .* newImagePatch));
        DiscStim = equivalentContrast .* centerBinary; %uniform disc
        
        %activation in space:
        convolved_SubunitCenter = conv2(DiscStim, SubunitFilter, 'same');
        subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
        convolved_SubunitSurround = conv2(DiscStim, SubunitSurroundFilter, 'same');
        subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
        %temporal activation of each subunit center & surround
        %rows = subunit, cols = time
        subunitVoltage = subunitCenterActivations * centerFilter -...
            subunitSuroundActivations * surroundFilter;
        %hit each subunit with rectifying synapse
        rectInds = find(subunitVoltage < 0);
        subunitVoltage(rectInds) = 0;
        subunitOutputs = subunitVoltage;
        response.Disc(pp,:) = sum(subunitOutputs .* subunitWeightings);


    %3) ImageMatchedSurround
        ImageMatchedSurroundStim = newImagePatch;
        %activation in space:
        convolved_SubunitCenter = conv2(ImageMatchedSurroundStim, SubunitFilter, 'same');
        subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
        convolved_SubunitSurround = conv2(ImageMatchedSurroundStim, SubunitSurroundFilter, 'same');
        subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
        %temporal activation of each subunit center & surround
        %rows = subunit, cols = time
        subunitVoltage = subunitCenterActivations * centerFilter -...
            subunitSuroundActivations * surroundFilter;
        %hit each subunit with rectifying synapse
        rectInds = find(subunitVoltage < 0);
        subunitVoltage(rectInds) = 0;
        subunitOutputs = subunitVoltage;
        response.ImageMatchedSurround(pp,:) = sum(subunitOutputs .* subunitWeightings);

    %4) DiscMatchedSurround
        DiscMatchedSurroundStim = newImagePatch;
        DiscMatchedSurroundStim(logical(centerBinary)) = equivalentContrast; %put equiv. disc in center
        %activation in space:
        convolved_SubunitCenter = conv2(DiscMatchedSurroundStim, SubunitFilter, 'same');
        subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
        convolved_SubunitSurround = conv2(DiscMatchedSurroundStim, SubunitSurroundFilter, 'same');
        subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
        %temporal activation of each subunit center & surround
        %rows = subunit, cols = time
        subunitVoltage = subunitCenterActivations * centerFilter -...
            subunitSuroundActivations * surroundFilter;
        %hit each subunit with rectifying synapse
        rectInds = find(subunitVoltage < 0);
        subunitVoltage(rectInds) = 0;
        subunitOutputs = subunitVoltage;
        response.DiscMatchedSurround(pp,:) = sum(subunitOutputs .* subunitWeightings);

    
    %5) ImageMixedSurround
        ImageMixedSurroundStim = mixedImagePatch;
        ImageMixedSurroundStim(logical(centerBinary)) = ImageStim(logical(centerBinary));
        %activation in space:
        convolved_SubunitCenter = conv2(ImageMixedSurroundStim, SubunitFilter, 'same');
        subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
        convolved_SubunitSurround = conv2(ImageMixedSurroundStim, SubunitSurroundFilter, 'same');
        subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
        %temporal activation of each subunit center & surround
        %rows = subunit, cols = time
        subunitVoltage = subunitCenterActivations * centerFilter -...
            subunitSuroundActivations * surroundFilter;
        %hit each subunit with rectifying synapse
        rectInds = find(subunitVoltage < 0);
        subunitVoltage(rectInds) = 0;
        subunitOutputs = subunitVoltage;
        response.ImageMixedSurround(pp,:) = sum(subunitOutputs .* subunitWeightings);
    
    %6) DiscMixedSurround
        DiscMixedSurroundStim = mixedImagePatch;
        DiscMixedSurroundStim(logical(centerBinary)) = equivalentContrast;
        %activation in space:
        convolved_SubunitCenter = conv2(DiscMixedSurroundStim, SubunitFilter, 'same');
        subunitCenterActivations = convolved_SubunitCenter(SubunitIndices);
    
        convolved_SubunitSurround = conv2(DiscMixedSurroundStim, SubunitSurroundFilter, 'same');
        subunitSuroundActivations = convolved_SubunitSurround(SubunitIndices);
    
        %temporal activation of each subunit center & surround
        %rows = subunit, cols = time
        subunitVoltage = subunitCenterActivations * centerFilter -...
            subunitSuroundActivations * surroundFilter;
        %hit each subunit with rectifying synapse
        rectInds = find(subunitVoltage < 0);
        subunitVoltage(rectInds) = 0;
        subunitOutputs = subunitVoltage;
        response.DiscMixedSurround(pp,:) = sum(subunitOutputs .* subunitWeightings);
    
end


%%
figure(4); clf;
for pp = 1:50
    subplot(211)
    plot(response.Image(pp,:),'k')
    hold on;
    plot(response.Disc(pp,:),'b')
    hold off;
    
    subplot(212)
    plot(response.ImageMatchedSurround(pp,:),'k')
    hold on;
    plot(response.DiscMatchedSurround(pp,:),'b')
    hold off;
    
    pause;
end

%%

figure(5); clf;
subplot(221);
intResponse_image = trapz(response.Image,2);
intResponse_disc = trapz(response.Disc,2);
respDiff_none = intResponse_image - intResponse_disc;
plot(intResponse_image,intResponse_disc,'ko')
hold on; plot([0, 300],[0, 300],'k--')

subplot(222);
intResponse_image = trapz(response.ImageMatchedSurround,2);
intResponse_disc = trapz(response.DiscMatchedSurround,2);
respDiff_nat = intResponse_image - intResponse_disc;
plot(intResponse_image,intResponse_disc,'go')
hold on; plot([0, 300],[0, 300],'k--')

subplot(223);
intResponse_image = trapz(response.ImageMixedSurround,2);
intResponse_disc = trapz(response.DiscMixedSurround,2);
respDiff_mix = intResponse_image - intResponse_disc;
plot(intResponse_image,intResponse_disc,'go')
hold on; plot([0, 300],[0, 300],'k--')

figure(6); clf;
subplot(221); plot(respDiff_none,respDiff_nat,'go')
hold on; plot([0, 100],[0, 100],'k--')

subplot(222); plot(respDiff_none,respDiff_mix,'ro')
hold on; plot([0, 100],[0, 100],'k--')

subplot(223); plot(respDiff_nat,respDiff_mix,'ko')
hold on; plot([0, 100],[0, 100],'k--')








