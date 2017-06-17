clear all; close all; clc;

stimInd = 1; %image to use
noPatches = 5000;
micronsPerPixel = 6.6;

%RF components:
SubunitSurroundWeight = 0.72;     % relative to center integral

filterSize_um = 700;                % size of patch (um)
subunitSigma_um = 10;
subunitSurroundSigma_um = 150;
centerSigma_um = 40;

resources_dir = '~/Documents/MATLAB/turner-package/resources/';
IMAGES_DIR_VH = '~/Documents/MATLAB/MHT-analysis/resources/vanhateren_iml/';
load([resources_dir, 'dovesFEMstims_20160422.mat'])

%image:
fileName = [FEMdata(stimInd).ImageName(1:8), '_',num2str(FEMdata(stimInd).SubjectIndex),'.mat'];
f1=fopen([IMAGES_DIR_VH, FEMdata(stimInd).ImageName],'rb','ieee-be');
w=1536;h=1024;
my_image=fread(f1,[w,h],'uint16');
my_image = my_image';

my_image_nomean = (my_image - mean(my_image(:))) ./ mean(my_image(:));

%convert to pixels:
filterSize = round(filterSize_um / micronsPerPixel);                % size of patch (um)
subunitSigma = round(subunitSigma_um / micronsPerPixel);
subunitSurroundSigma = round(subunitSurroundSigma_um / micronsPerPixel);
centerSigma = round(centerSigma_um / micronsPerPixel);

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


figure(2); clf;  set(gcf, 'WindowStyle', 'docked')
subplot(221); imagesc(RFCenter); colormap(gray); colorbar; title('Center')
subplot(222); imagesc(SubunitWithSurroundFilter); colormap(gray); colorbar; title('Subunit w/ surround')

% response to spots

stimulus.spotSize = 0:10:700; %diameter of spot

response.spots = [];

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
for ss = 1:length(stimulus.spotSize) %get responses to each spot
    currentRadius = (stimulus.spotSize(ss)/2)/ micronsPerPixel; %convert to pixel
    spotBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=currentRadius);

    %Model: Subunits have surrounds:
    convolved_SubunitWithSurround = conv2(spotBinary, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
    subunitActivations = convolved_SubunitWithSurround(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.spots(ss) = sum(subunitOutputs .* subunitWeightings);
end

[~, ind] = max(response.spots);
centerSize_um = stimulus.spotSize(ind);
centerSize = round(centerSize_um / micronsPerPixel);


%Exp spots responses
figure; clf;
fig1=gca;
set(fig1,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Spot diameter (um)')
set(get(fig1,'YLabel'),'String','Response (a.u.)')
addLineToAxis(stimulus.spotSize,response.spots,'areaSum',fig1,'k','-','none')

figID = 'MixedSurModel_ExpSpots';
makeAxisStruct(fig1,figID ,'RFSurroundFigs')


%Subunit filter
figure; clf;
fig3=gca;
set(fig3,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig3,'XLabel'),'String','Loc (um)')
set(get(fig3,'YLabel'),'String','Sens')
addLineToAxis((-filterSize/2 + 1) : (filterSize/2),sum(SubunitWithSurroundFilter),'areaSum',fig3,'k','-','none')

figID = 'MixedSurModel_SubFilter';
makeAxisStruct(fig3,figID ,'RFSurroundFigs')
%% responses to image patches and discs, +/- surrounds
tic;
contrastPolarity = -1;

%center only:
response.Image = [];
response.Disc = [];
%matched surround:
response.ImageMatchedSurround = [];
response.DiscMatchedSurround = [];
%mixed surround:
response.ImageMixedSurround = [];
response.DiscMixedSurround = [];

stimulus.CenterLocation = [];
stimulus.MixedSurroundLocation = [];


[rr, cc] = meshgrid(1:filterSize,1:filterSize);
centerBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=(centerSize/2));
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
    
    %1) Image:
        ImageStim = newImagePatch;
        ImageStim(~centerBinary) = 0; %add aperture

        convolved_im = conv2(ImageStim, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
        subunitActivations = convolved_im(SubunitIndices);
        subunitOutputs = subunitActivations;
        subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        response.Image(pp) = sum(subunitOutputs .* subunitWeightings);

    %2) Disc:
        %compute linear equivalent contrast...
        % Get the model RF:
        RF = fspecial('gaussian',filterSize,centerSigma);
        weightingFxn = centerBinary .* RF; %set to zero mean gray pixels
        weightingFxn = weightingFxn ./ sum(weightingFxn(:)); %sum to one
        equivalentContrast = sum(sum(weightingFxn .* newImagePatch));
        DiscStim = equivalentContrast .* centerBinary; %uniform disc
        
        convolved_im = conv2(DiscStim, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
        subunitActivations = convolved_im(SubunitIndices);
        subunitOutputs = subunitActivations;
        subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        response.Disc(pp) = sum(subunitOutputs .* subunitWeightings);


    %3) ImageMatchedSurround
        ImageMatchedSurroundStim = newImagePatch;
        
        convolved_im = conv2(ImageMatchedSurroundStim, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
        subunitActivations = convolved_im(SubunitIndices);
        subunitOutputs = subunitActivations;
        subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        response.ImageMatchedSurround(pp) = sum(subunitOutputs .* subunitWeightings);
    
    %4) DiscMatchedSurround
        DiscMatchedSurroundStim = newImagePatch;
        DiscMatchedSurroundStim(logical(centerBinary)) = equivalentContrast; %put equiv. disc in center
        
        convolved_im = conv2(DiscMatchedSurroundStim, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
        subunitActivations = convolved_im(SubunitIndices);
        subunitOutputs = subunitActivations;
        subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        response.DiscMatchedSurround(pp) = sum(subunitOutputs .* subunitWeightings);
    
    
    %5) ImageMixedSurround
        ImageMixedSurroundStim = mixedImagePatch;
        ImageMixedSurroundStim(logical(centerBinary)) = ImageStim(logical(centerBinary));
        
        convolved_im = conv2(ImageMixedSurroundStim, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
        subunitActivations = convolved_im(SubunitIndices);
        subunitOutputs = subunitActivations;
        subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        response.ImageMixedSurround(pp) = sum(subunitOutputs .* subunitWeightings);
    
    %6) DiscMixedSurround
        DiscMixedSurroundStim = mixedImagePatch;
        DiscMixedSurroundStim(logical(centerBinary)) = equivalentContrast;
        
        convolved_im = conv2(DiscMixedSurroundStim, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
        subunitActivations = convolved_im(SubunitIndices);
        subunitOutputs = subunitActivations;
        subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        response.DiscMixedSurround(pp) = sum(subunitOutputs .* subunitWeightings);
    
%     figDir = '~/Documents/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs
%     fh = figure(30); clf;
%     upLim = max(ImageMatchedSurroundStim(:)); downLim = min(ImageMatchedSurroundStim(:));
%     subplot(321); imagesc(contrastPolarity .* ImageStim); colormap(gray); caxis([downLim, upLim]); axis image; axis off; 
%     subplot(322); imagesc(contrastPolarity .* DiscStim); colormap(gray); caxis([downLim, upLim]); axis image; axis off; 
%     
%     subplot(323); imagesc(contrastPolarity .* ImageMatchedSurroundStim); colormap(gray); caxis([downLim, upLim]); axis image; axis off; 
%     subplot(324); imagesc(contrastPolarity .* DiscMatchedSurroundStim); colormap(gray); caxis([downLim, upLim]); axis image; axis off; 
%     
%     subplot(325); imagesc(contrastPolarity .* ImageMixedSurroundStim); colormap(gray); caxis([downLim, upLim]); axis image; axis off; 
%     subplot(326); imagesc(contrastPolarity .* DiscMixedSurroundStim); colormap(gray); caxis([downLim, upLim]); axis image; axis off; 
% %     brighten(0.6) %brighten for display purposes
%     drawnow;
%     figID = 'egStims_MixedSurround';
%     print(fh,[figDir,figID],'-depsc')

    
   
end
toc;

save(['ImageDiscModelResults_',num2str(stimInd),'_20170616.mat'],'stimulus','response')

%%


NLI_noSurround = (response.Image - response.Disc) ./...
    (response.Image + response.Disc);

NLI_naturalSurround = (response.ImageMatchedSurround - response.DiscMatchedSurround) ./...
    (response.ImageMatchedSurround + response.DiscMatchedSurround);

NLI_mixedSurround = (response.ImageMixedSurround - response.DiscMixedSurround) ./...
    (response.ImageMixedSurround + response.DiscMixedSurround);

[val, ind] = sort(NLI_noSurround, 'descend');

figure(20); clf; subplot(211);
hold on;
plot(1:noPatches,val,'k.')
plot(1:noPatches,NLI_naturalSurround(ind),'ro')

subplot(212);
hold on;
plot(1:noPatches,val,'k.')
plot(1:noPatches,NLI_mixedSurround(ind),'ro')

cutoffNLI = 0.1;
takeUpTo = find(val < cutoffNLI,1) - 1; %cutoff at cutoffNLI
noBins = 25;

%NLI histograms
figure; clf;
fig2=gca;
set(fig2,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig2,'XLabel'),'String','NLI')
set(get(fig2,'YLabel'),'String','Fraction image patches')
addLineToAxis([cutoffNLI cutoffNLI],[0 0.4],'cutoff',fig2,'k','--','none')
%no surround:
[nn, ctr] = histcounts(val(1:takeUpTo),noBins);
binCtrs = ctr(1:end-1) + diff(ctr);
addLineToAxis(binCtrs,nn./sum(nn),'noSurround',fig2,'k','-','none')
%nat surround:
[nn, ctr] = histcounts(NLI_naturalSurround(ind(1:takeUpTo)),noBins);
binCtrs = ctr(1:end-1) + diff(ctr);
addLineToAxis(binCtrs,nn./sum(nn),'natSurround',fig2,'g','-','none')
%mixed surround:
[nn, ctr] = histcounts(NLI_mixedSurround(ind(1:takeUpTo)),noBins);
binCtrs = ctr(1:end-1) + diff(ctr);
addLineToAxis(binCtrs,nn./sum(nn),'mixSurround',fig2,'r','-','none')

figID = 'MixedSurModel_NLIhist';
makeAxisStruct(fig2,figID ,'RFSurroundFigs')


