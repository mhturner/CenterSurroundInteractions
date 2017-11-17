%% % % % % % % % % RESPONSE TO SPOTS % % % % % % % % % % % % % % % % % % % % % 
clear all; close all; clc;
RF = SpatialRFModel;
RF.subunitSurroundWeight = 0.4; %0.72 ~Off parasol; 0.9 = str; 0.4 = wk
RF.makeRfComponents;

stimulus.spotSize = 0:10:700; %diameter of spot

response.spots = [];
filterSize = size(RF.SubunitFilter,1);

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
colors = pmkmp(length(stimulus.spotSize));
for ss = 1:length(stimulus.spotSize) %get responses to each spot
    currentRadius = (stimulus.spotSize(ss)/2)/ RF.MicronsPerPixel; %convert to pixel
    spotBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=currentRadius);
    
    response.spots(ss) = RF.getResponse(spotBinary);
end

figure; clf; fig1=gca; initFig(fig1,'Spot size (um)','Resp. (au)')
addLineToAxis(stimulus.spotSize,response.spots,'dat',fig1,'k','none','o')

[~, ind] = max(response.spots);
centerSize_um = stimulus.spotSize(ind);
centerSize = round(centerSize_um / RF.MicronsPerPixel);
makeAxisStruct(fig1,'SModel_AS_wk' ,'RFSurroundFigs')

figure; clf; fig2=gca; initFig(fig2,'loc','Resp. (au)')
xx = 1:size(RF.SubunitFilter,1);
addLineToAxis(xx,mean(RF.SubunitFilter),'Cfilt',fig2,'g','-','none')
addLineToAxis(xx,-mean(RF.SubunitSurroundFilter),'Sfilt',fig2,'m','-','none')
addLineToAxis(xx,mean(RF.SubunitFilter) - mean(RF.SubunitSurroundFilter),'Sumfilt',fig2,'k','-','none')
makeAxisStruct(fig2,'SModel_SubFilters' ,'RFSurroundFigs')

%% % % % % % % % % RESPONSE TO NATURAL IMAGE STIMS % % % % % % % % % % % % % % % % % % % %
noPatches = 50;
noImages = 29;

RF = SpatialRFModel;
RF.subunitSurroundWeight = 0.4; %0.72 ~Off parasol; 0.9 = str; 0.4 = wk
RF.makeRfComponents;

filterSize = size(RF.SubunitFilter,1);
[rr, cc] = meshgrid(1:filterSize,1:filterSize);
centerBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=(centerSize/2));
surroundBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)>(centerSize/2));

IMAGES_DIR_VH = '~/Dropbox/RiekeLab/Analysis/MATLAB/MHT-analysis/resources/vanhateren_iml/';
resources_dir = '~/Dropbox/RiekeLab/Analysis/MATLAB/turner-package/resources/';
load([resources_dir, 'dovesFEMstims_20160422.mat'])

allNLI.none = [];
allNLI.nat = [];
allNLI.mix = [];
allPatchIntensity = [];
for ii = 1:noImages
targetImageIndex = ii;
% % % % % % % % LOAD NATURAL IMAGE % % % %
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

patchIntensity = [];
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
        response.Image(pp,:) = RF.getResponse(ImageStim);
    %2) Disc:
        %compute linear equivalent contrast...
        % Get the model RF:
        centerRF = fspecial('gaussian',filterSize,(RF.centerSigma_um / RF.MicronsPerPixel));
        weightingFxn = centerBinary .* centerRF; %set to zero mean gray pixels
        weightingFxn = weightingFxn ./ sum(weightingFxn(:)); %sum to one
        equivalentContrast = sum(sum(weightingFxn .* newImagePatch));
        DiscStim = equivalentContrast .* centerBinary; %uniform disc
        response.Disc(pp,:) = RF.getResponse(DiscStim);
    %3) ImageMatchedSurround
        ImageMatchedSurroundStim = newImagePatch;
        response.ImageMatchedSurround(pp,:) = RF.getResponse(ImageMatchedSurroundStim);
    %4) DiscMatchedSurround
        DiscMatchedSurroundStim = newImagePatch;
        DiscMatchedSurroundStim(logical(centerBinary)) = equivalentContrast; %put equiv. disc in center
        response.DiscMatchedSurround(pp,:) = RF.getResponse(DiscMatchedSurroundStim);
    %5) ImageMixedSurround
        ImageMixedSurroundStim = mixedImagePatch;
        ImageMixedSurroundStim(logical(centerBinary)) = ImageStim(logical(centerBinary));
        response.ImageMixedSurround(pp,:) = RF.getResponse(ImageMixedSurroundStim);
    %6) DiscMixedSurround
        DiscMixedSurroundStim = mixedImagePatch;
        DiscMixedSurroundStim(logical(centerBinary)) = equivalentContrast;
        response.DiscMixedSurround(pp,:) = RF.getResponse(DiscMixedSurroundStim);
        
        %eq contrast calc based on contrast polarity flipped. Flip back...
        patchIntensity(pp) = contrastPolarity .* equivalentContrast;
end

allNLI.none(ii,:) = (response.Image - response.Disc) ./ ...
    (response.Image + response.Disc);

allNLI.nat(ii,:) = (response.ImageMatchedSurround - response.DiscMatchedSurround) ./ ...
    (response.ImageMatchedSurround + response.DiscMatchedSurround);

allNLI.mix(ii,:) = (response.ImageMixedSurround - response.DiscMixedSurround) ./ ...
    (response.ImageMixedSurround + response.DiscMixedSurround);

allPatchIntensity(ii,:) = patchIntensity;
disp(ii)
end
%%
NLInone = allNLI.none(:);
NLInat = allNLI.nat(:);
NLImix = allNLI.mix(:);
patchIntensity = allPatchIntensity(:);
cutInds = find(patchIntensity > 0); %cut positive means

NLInone(cutInds) = [];
NLInat(cutInds) = [];
NLImix(cutInds) = [];
patchIntensity(cutInds) = [];

noBins = 12;

figure; clf; fig11=gca; initFig(fig11,'Relative center mean','NLI') % NLI vs center intensity

tempNLIs = NLInone;
binAndPlotEquallyPopulatedBins(patchIntensity(~isnan(tempNLIs)),tempNLIs(~isnan(tempNLIs)),noBins,fig11,'k','none')

tempNLIs = NLInat;
binAndPlotEquallyPopulatedBins(patchIntensity(~isnan(tempNLIs)),tempNLIs(~isnan(tempNLIs)),noBins,fig11,'g','nat')

makeAxisStruct(fig11,'SModel_NLIbinned_wk' ,'RFSurroundFigs')

