%% % % % % % % % % RESPONSE TO SPOTS % % % % % % % % % % % % % % % % % % % % % 
clear all; close all; clc;
RF = SpatioTemporalRFModel;
RF.subunitSurroundWeight = 0.5; %1, at 0 shift ~Off parasol; 1.5 = str; 0.5 = wk
RF.surroundFilterTimeShift = 0;
stimulusStep = [zeros(1,150),ones(1,200),zeros(1,150)];
timeVec = (0:499);

RF.makeRfComponents;

stimulus.spotSize = 0:10:700; %diameter of spot

response.spots = [];
filterSize = size(RF.SubunitFilter,1);

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
colors = pmkmp(length(stimulus.spotSize));

figure; clf; fig1=gca; initFig(fig1,'Time (ms)','Resp. (au)')
for ss = 1:length(stimulus.spotSize) %get responses to each spot
    currentRadius = (stimulus.spotSize(ss)/2)/ RF.MicronsPerPixel; %convert to pixel
    spotBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=currentRadius);
    
    response.spots(ss,:) = RF.getResponse(spotBinary,stimulusStep);
    addLineToAxis(timeVec,response.spots(ss,:),['resp',num2str(ss)],fig1,colors(ss,:),'-','none')
end

intResponse = trapz(response.spots(:,150:350),2);
figure; clf; fig2=gca; initFig(fig2,'Spot size (um)','Resp. (au)')
addLineToAxis(stimulus.spotSize,intResponse,'dat',fig2,'k','none','o')


figure; clf; fig3=gca; initFig(fig3,'Time (ms)','Resp. (au)')
addLineToAxis(timeVec,RF.CenterTemporalFilter,'Cfilt',fig3,'g','-','none')
addLineToAxis(timeVec,RF.SurroundTemporalFilter,'Sfilt',fig3,'m','-','none')


[~, ind] = max(intResponse);
centerSize_um = stimulus.spotSize(ind);
centerSize = round(centerSize_um / RF.MicronsPerPixel);

makeAxisStruct(fig1,'TModel_ES' ,'RFSurroundFigs')
makeAxisStruct(fig2,'TModel_AS_wk' ,'RFSurroundFigs')
makeAxisStruct(fig3,'TModel_filters' ,'RFSurroundFigs')

figure; clf; fig4=gca; initFig(fig4,'loc','Resp. (au)')
xx = 1:size(RF.SubunitFilter,1);
addLineToAxis(xx,mean(RF.SubunitFilter),'Cfilt',fig4,'g','-','none')
addLineToAxis(xx,-mean(RF.SubunitSurroundFilter),'Sfilt',fig4,'m','-','none')
addLineToAxis(xx,mean(RF.SubunitFilter) - mean(RF.SubunitSurroundFilter),'Sumfilt',fig4,'k','-','none')
makeAxisStruct(fig4,'TModel_SubFilters' ,'RFSurroundFigs')

%% % % % % % % % % RESPONSE TO NATURAL IMAGE STIMS % % % % % % % % % % % % % % % % % % % %
noPatches = 50;
noImages = 29;


RF = SpatioTemporalRFModel;
RF.subunitSurroundWeight = 0.5; %1, at 0 shift ~Off parasol; 1.5 = str; 0.5 = wk;
RF.surroundFilterTimeShift = 0; %negative means move surround left (faster)
stimulusStep = [zeros(1,150),ones(1,200),zeros(1,150)];
intInds = 150:350;

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
        response.Image(pp,:) = RF.getResponse(ImageStim,stimulusStep);
    %2) Disc:
        %compute linear equivalent contrast...
        % Get the model RF:
        centerRF = fspecial('gaussian',filterSize,(RF.centerSigma_um / RF.MicronsPerPixel));
        weightingFxn = centerBinary .* centerRF; %set to zero mean gray pixels
        weightingFxn = weightingFxn ./ sum(weightingFxn(:)); %sum to one
        equivalentContrast = sum(sum(weightingFxn .* newImagePatch));
        DiscStim = equivalentContrast .* centerBinary; %uniform disc
        response.Disc(pp,:) = RF.getResponse(DiscStim,stimulusStep);
    %3) ImageMatchedSurround
        ImageMatchedSurroundStim = newImagePatch;
        response.ImageMatchedSurround(pp,:) = RF.getResponse(ImageMatchedSurroundStim,stimulusStep);
    %4) DiscMatchedSurround
        DiscMatchedSurroundStim = newImagePatch;
        DiscMatchedSurroundStim(logical(centerBinary)) = equivalentContrast; %put equiv. disc in center
        response.DiscMatchedSurround(pp,:) = RF.getResponse(DiscMatchedSurroundStim,stimulusStep);
    %5) ImageMixedSurround
        ImageMixedSurroundStim = mixedImagePatch;
        ImageMixedSurroundStim(logical(centerBinary)) = ImageStim(logical(centerBinary));
        response.ImageMixedSurround(pp,:) = RF.getResponse(ImageMixedSurroundStim,stimulusStep);
    %6) DiscMixedSurround
        DiscMixedSurroundStim = mixedImagePatch;
        DiscMixedSurroundStim(logical(centerBinary)) = equivalentContrast;
        response.DiscMixedSurround(pp,:) = RF.getResponse(DiscMixedSurroundStim,stimulusStep);
        
        %eq contrast calc based on contrast polarity flipped. Flip back...
        patchIntensity(pp) = contrastPolarity .* equivalentContrast;
end

intResponse_image = trapz(response.Image(:,intInds),2);
intResponse_disc = trapz(response.Disc(:,intInds),2);
allNLI.none(ii,:) = (intResponse_image - intResponse_disc) ./ ...
    (intResponse_image + intResponse_disc);

intResponse_image = trapz(response.ImageMatchedSurround(:,intInds),2);
intResponse_disc = trapz(response.DiscMatchedSurround(:,intInds),2);
allNLI.nat(ii,:) = (intResponse_image - intResponse_disc) ./ ...
    (intResponse_image + intResponse_disc);

intResponse_image = trapz(response.ImageMixedSurround(:,intInds),2);
intResponse_disc = trapz(response.DiscMixedSurround(:,intInds),2);
allNLI.mix(ii,:) = (intResponse_image - intResponse_disc) ./ ...
    (intResponse_image + intResponse_disc);

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

makeAxisStruct(fig11,'TModel_NLIbinned_wk' ,'RFSurroundFigs')

