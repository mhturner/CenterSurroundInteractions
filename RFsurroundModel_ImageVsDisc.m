%% Compute new RF surround model results
clear all; close all; clc;

outputDir = '~/Documents/MATLAB/RFSurround/resources/RFsurroundNLIresults/';
resources_dir = '~/Documents/MATLAB/turner-package/resources/';
IMAGES_DIR_VH = '~/Documents/MATLAB/MHT-analysis/resources/vanhateren_iml/';
load([resources_dir, 'dovesFEMstims_20160422.mat'])

for currentImageIndex = 1:30

targetImageIndex = currentImageIndex; %up to 30
for ss = 1:length(FEMdata)
    if FEMdata(ss).ImageIndex == targetImageIndex;
        stimInd = ss;
        break;
    end
end

noPatches = 1000;
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

% 
% figure(2); clf;  set(gcf, 'WindowStyle', 'docked')
% subplot(221); imagesc(RFCenter); colormap(gray); colorbar; title('Center')
% subplot(222); imagesc(SubunitWithSurroundFilter); colormap(gray); colorbar; title('Subunit w/ surround')

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
% figure; clf;
% fig1=gca;
% set(fig1,'XScale','linear','YScale','linear')
% set(0, 'DefaultAxesFontSize', 12)
% set(get(fig1,'XLabel'),'String','Spot diameter (um)')
% set(get(fig1,'YLabel'),'String','Response (a.u.)')
% addLineToAxis(stimulus.spotSize,response.spots,'areaSum',fig1,'k','-','none')
% 
% figID = 'MixedSur_ES_wk';
% makeAxisStruct(fig1,figID ,'RFSurroundFigs')
% 
% 
% %Subunit filter
% figure; clf;
% fig3=gca;
% set(fig3,'XScale','linear','YScale','linear')
% set(0, 'DefaultAxesFontSize', 12)
% set(get(fig3,'XLabel'),'String','Loc (um)')
% set(get(fig3,'YLabel'),'String','Sens')
% addLineToAxis((-filterSize/2 + 1) : (filterSize/2),sum(SubunitWithSurroundFilter),'areaSum',fig3,'k','-','none')
% 
% figID = 'MixedSur_Sub_wk';
% makeAxisStruct(fig3,figID ,'RFSurroundFigs')
% % % responses to image patches and discs, +/- surrounds

%
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

stimulus.CenterLocation = [];
stimulus.MixedSurroundLocation = [];

stats.centerMean = [];
stats.centerVar = [];

stats.surroundMean = [];
stats.surroundVar = [];

stats.surroundMean_Mix = [];
stats.surroundVar_Mix = [];

stats.imageMean = mean(my_image(:));

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
centerBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=(centerSize/2));
surroundBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)>(centerSize/2));
for pp = 1:noPatches
    %main image patch location:
    x = round(filterSize/2 + (w - filterSize)*rand);
    y = round(filterSize/2 + (h - filterSize)*rand);
    newImagePatch = contrastPolarity .* my_image_nomean(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    stimulus.CenterLocation(pp,:) = [x, y];
    
    %get main location stats:
    statsImagePatch = my_image(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    
    centerPixels = statsImagePatch(find(centerBinary));
    stats.centerMean(pp) = mean(centerPixels);
    stats.centerVar(pp) = var(centerPixels);
    
    surroundPixels = statsImagePatch(find(surroundBinary));
    stats.surroundMean(pp) = mean(surroundPixels);
    stats.surroundVar(pp) = var(surroundPixels);
    
    %mixed image surround location:
    x_mix = round(filterSize/2 + (w - filterSize)*rand);
    y_mix = round(filterSize/2 + (h - filterSize)*rand);
    mixedImagePatch = contrastPolarity .* my_image_nomean(y_mix - filterSize/2 + 1 : y_mix + filterSize/2, x_mix - filterSize/2 + 1 : x_mix + filterSize/2);
    stimulus.MixedSurroundLocation(pp,:) = [x_mix, y_mix];
    
    %get mixed location stats:
    statsImagePatch = my_image(y_mix - filterSize/2 + 1 : y_mix + filterSize/2, x_mix - filterSize/2 + 1 : x_mix + filterSize/2);
    surroundPixels = statsImagePatch(find(surroundBinary));
    stats.surroundMean_Mix(pp) = mean(surroundPixels);
    stats.surroundVar_Mix(pp) = var(surroundPixels);

    %Stim responses:
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
% save([outputDir,'ImageDiscModel_',num2str(targetImageIndex),'_20170718.mat'],'stimulus','response','stats')

disp(currentImageIndex)

end

%% Center-surround intensity difference histogram and image
clear all; close all; clc;
egIm = 1;


resources_dir = '~/Documents/MATLAB/turner-package/resources/';
IMAGES_DIR_VH = '~/Documents/MATLAB/MHT-analysis/resources/vanhateren_iml/';
load([resources_dir, 'dovesFEMstims_20160422.mat'])
figDir = '~/Documents/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs

ct = 0;
results.natDiff = [];
results.mixDiff = [];


noBins = 40;
binEdges = linspace(-1,1,noBins+1);
binCtrs = binEdges(1:end-1) + diff(binEdges);

for currentImageIndex = 1:30
    
    if currentImageIndex == 12
        continue
    end
    load(['ImageDiscModel_',num2str(currentImageIndex),'_20170718.mat'])
ct = ct + 1;

intDiff_nat = (stats.centerMean - stats.surroundMean) ./ stats.imageMean;
intDiff_mix = (stats.centerMean - stats.surroundMean_Mix) ./ stats.imageMean;

[nn_nat, ~] = histcounts(intDiff_nat,binEdges,'Normalization','probability');
[nn_mix, ~] = histcounts(intDiff_mix,binEdges,'Normalization','probability');


    if currentImageIndex == egIm
        for ss = 1:length(FEMdata)
            if FEMdata(ss).ImageIndex == egIm;
                stimInd = ss;
                break;
            end
        end
        
        %image:
        fileName = [FEMdata(stimInd).ImageName(1:8), '_',num2str(FEMdata(stimInd).SubjectIndex),'.mat'];
        f1=fopen([IMAGES_DIR_VH, FEMdata(stimInd).ImageName],'rb','ieee-be');
        w=1536;h=1024;
        img=fread(f1,[w,h],'uint16');
        img = img';
        img = img ./ max(img(:)); %normalize to [0 1]
        fh = figure(21); clf;
        imagesc(img); colormap(gray); axis image; axis off; hold on;
        brighten(0.6) %brighten for display purposes
        drawnow;
        figID = 'egNatImage_MixSur';
        print(fh,[figDir,figID],'-depsc')
        

        figure; clf;
        fig4=gca;
        set(fig4,'XScale','linear','YScale','linear')
        set(0, 'DefaultAxesFontSize', 12)
        set(get(fig4,'XLabel'),'String','Ic - Is')
        set(get(fig4,'YLabel'),'String','Prob.')
        addLineToAxis(binCtrs,nn_nat,'nat',fig4,'k','-','none')
        addLineToAxis(binCtrs,nn_mix,'mix',fig4,'r','-','none')
    end
    results.natDiff = cat(1,results.natDiff,nn_nat);
    results.mixDiff = cat(1,results.mixDiff,nn_mix);

end
figure; clf;
fig5=gca;
set(fig5,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig5,'XLabel'),'String','Ic - Is')
set(get(fig5,'YLabel'),'String','Prob.')

errNat = std(results.natDiff,[],1) ./ sqrt(size(results.natDiff,1));
addLineToAxis(binCtrs,mean(results.natDiff,1),'natMean',fig5,'k','-','none')
addLineToAxis(binCtrs,mean(results.natDiff,1) + errNat,'natErrUp',fig5,'k','--','none')
addLineToAxis(binCtrs,mean(results.natDiff,1) - errNat,'natErrDown',fig5,'k','--','none')

errMix = std(results.mixDiff,[],1) ./ sqrt(size(results.mixDiff,1));
addLineToAxis(binCtrs,mean(results.mixDiff,1),'mixMean',fig5,'r','-','none')
addLineToAxis(binCtrs,mean(results.mixDiff,1) + errMix,'mixErrUp',fig5,'r','--','none')
addLineToAxis(binCtrs,mean(results.mixDiff,1) - errMix,'mixErrDown',fig5,'r','--','none')

figID = ['MixSurStats_eg'];
makeAxisStruct(fig4,figID ,'RFSurroundFigs')

figID = ['MixSurStats_pop'];
makeAxisStruct(fig5,figID ,'RFSurroundFigs')


%% R-squared between image and disc responses under three surround conditions:

noPts = 50;

figure(15); clf;
ct = 0;
p_none = [];
p_nat = [];
p_mix = [];
for currentImageIndex = 1:30
    
    if currentImageIndex == 12
        continue
    end
    ct = ct + 1;
    load(['ImageDiscModel_',num2str(currentImageIndex),'_20170718.mat'])

    
% %     ssTot=sum((response.Image-mean(response.Image)).^2); 
% %     ssErr=sum((response.Image - response.Disc).^2);
% %     rSquared=1-ssErr/ssTot;
% %     p_none(ct) = rSquared;
% %     
% %     ssTot=sum((response.ImageMatchedSurround-mean(response.ImageMatchedSurround)).^2); 
% %     ssErr=sum((response.ImageMatchedSurround - response.DiscMatchedSurround).^2);
% %     rSquared=1-ssErr/ssTot;
% %     p_nat(ct) = rSquared;
% %     
% %     ssTot=sum((response.ImageMixedSurround-mean(response.ImageMixedSurround)).^2); 
% %     ssErr=sum((response.ImageMixedSurround - response.DiscMixedSurround).^2);
% %     rSquared=1-ssErr/ssTot;
% %     p_mix(ct) = rSquared;

    pickInds = randsample(1:length(response.Image),noPts);
    response.Image = response.Image(pickInds);
    response.Disc = response.Disc(pickInds);
    ssTot=sum((response.Image-mean(response.Image)).^2); 
    ssErr=sum((response.Image - response.Disc).^2);
    rSquared=1-ssErr/ssTot;
    p_none(ct) = rSquared;
    
    response.ImageMatchedSurround = response.ImageMatchedSurround(pickInds);
    response.DiscMatchedSurround = response.DiscMatchedSurround(pickInds);
    ssTot=sum((response.ImageMatchedSurround-mean(response.ImageMatchedSurround)).^2); 
    ssErr=sum((response.ImageMatchedSurround - response.DiscMatchedSurround).^2);
    rSquared=1-ssErr/ssTot;
    p_nat(ct) = rSquared;
    
    response.ImageMixedSurround = response.ImageMixedSurround(pickInds);
    response.DiscMixedSurround = response.DiscMixedSurround(pickInds);
    ssTot=sum((response.ImageMixedSurround-mean(response.ImageMixedSurround)).^2); 
    ssErr=sum((response.ImageMixedSurround - response.DiscMixedSurround).^2);
    rSquared=1-ssErr/ssTot;
    p_mix(ct) = rSquared;

end

subplot(131); hold on;
plot(p_none,p_nat,'go')
plot([0 1],[0 1],'k--')
xlabel('R2, no surround'); ylabel('R2, nat surround');
xlim([0 1]); ylim([0 1]);

subplot(132); hold on;
plot(p_none,p_mix,'ro')
plot([0 1],[0 1],'k--')
xlabel('R2, no surround'); ylabel('R2, mixed surround');
% xlim([0 1]); ylim([0 1]);

subplot(133); hold on;
plot(p_nat,p_mix,'ko')
plot([0 1],[0 1],'k--')
xlabel('R2, nat surround'); ylabel('R2, mixed surround');
% xlim([0 1]); ylim([0 1]);


[h, p] = ttest(p_none,p_nat);
[h, p] = ttest(p_none,p_mix);
[h, p] = ttest(p_nat,p_mix);


%% mean NLIs across all patches within an image, for different surround conditions
clear all; close all; clc;

surroundType = 'std';
% surroundType = 'weak';
% surroundType = 'strong';

noPatches = 1000;
noBins = 100;
binEdges = linspace(-1,1,noBins+1);
binCtrs = binEdges(1:end-1) + diff(binEdges);

all.nn.none = [];
all.nn.nat = [];
all.nn.mix = [];

low.nn.none = [];
low.nn.nat = [];
low.nn.mix = [];

high.nn.none = [];
high.nn.nat = [];
high.nn.mix = [];

meanNLI_none = [];
meanNLI_nat = [];
meanNLI_mix = [];

%skip 12, no image
for currentImageIndex = 1:30
    if currentImageIndex == 12
        continue
    end
    if strcmp(surroundType,'std')
        load(['ImageDiscModelResults_',num2str(currentImageIndex),'_20170623.mat'])
        figTag = '';
    elseif strcmp(surroundType,'weak')
        load(['ImageDiscModel_weak_',num2str(currentImageIndex),'_20170705.mat'])
        figTag = 'wk';
    elseif strcmp(surroundType,'strong')
        load(['ImageDiscModel_strong_',num2str(currentImageIndex),'_20170705.mat'])
        figTag = 'str';
    end
    

    %calc NLIs:
    NLI_noSurround = (response.Image - response.Disc) ./...
        (response.Image + response.Disc);
    
    NLI_naturalSurround = (response.ImageMatchedSurround - response.DiscMatchedSurround) ./...
        (response.ImageMatchedSurround + response.DiscMatchedSurround);

    NLI_mixedSurround = (response.ImageMixedSurround - response.DiscMixedSurround) ./...
        (response.ImageMixedSurround + response.DiscMixedSurround);
    
    meanNLI_none = cat(1,meanNLI_none,nanmean(NLI_noSurround));
    meanNLI_nat = cat(1,meanNLI_nat,nanmean(NLI_naturalSurround));
    meanNLI_mix = cat(1,meanNLI_mix,nanmean(NLI_mixedSurround));
    
    cutoffValue = mean(NLI_noSurround);
    
    %NLI based cutoff, both sides:
    lowIndices = find(NLI_noSurround < cutoffValue);
    highIndices = find(NLI_noSurround >= cutoffValue);
    
    %histograms:
    %all:
    [tempnn, ~] = histcounts(NLI_noSurround,binEdges,'Normalization','probability');
    all.nn.none = cat(1,all.nn.none,tempnn);
    
    [tempnn, ~] = histcounts(NLI_naturalSurround,binEdges,'Normalization','probability');
    all.nn.nat = cat(1,all.nn.nat,tempnn);
    
    [tempnn, ~] = histcounts(NLI_mixedSurround,binEdges,'Normalization','probability');
    all.nn.mix = cat(1,all.nn.mix,tempnn);
    
    
    %low starting NLI:
    [tempnn, ~] = histcounts(NLI_noSurround(lowIndices),binEdges,'Normalization','probability');
    low.nn.none = cat(1,low.nn.none,tempnn);
    
    [tempnn, ~] = histcounts(NLI_naturalSurround(lowIndices),binEdges,'Normalization','probability');
    low.nn.nat = cat(1,low.nn.nat,tempnn);
    
    [tempnn, ~] = histcounts(NLI_mixedSurround(lowIndices),binEdges,'Normalization','probability');
    low.nn.mix = cat(1,low.nn.mix,tempnn);
    
    %high starting NLI:
    [tempnn, ~] = histcounts(NLI_noSurround(highIndices),binEdges,'Normalization','probability');
    high.nn.none = cat(1,high.nn.none,tempnn);
    
    [tempnn, ~] = histcounts(NLI_naturalSurround(highIndices),binEdges,'Normalization','probability');
    high.nn.nat = cat(1,high.nn.nat,tempnn);
    
    [tempnn, ~] = histcounts(NLI_mixedSurround(highIndices),binEdges,'Normalization','probability');
    high.nn.mix = cat(1,high.nn.mix,tempnn);

end

figure; clf;
fig4=gca;
set(fig4,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig4,'XLabel'),'String','NLI no surround')
set(get(fig4,'YLabel'),'String','NLI natural surround')
addLineToAxis(meanNLI_none,meanNLI_nat,'data',fig4,'g','none','o')
errY = std(meanNLI_nat) ./ sqrt(length(meanNLI_nat));
errX = std(meanNLI_none) ./ sqrt(length(meanNLI_none));

addLineToAxis(mean(meanNLI_none),mean(meanNLI_nat),'mean',fig4,'g','none','s')

addLineToAxis([mean(meanNLI_none), mean(meanNLI_none)],...
    mean(meanNLI_nat) + [-errY, errY],'errY',fig4,'g','-','none')
addLineToAxis(mean(meanNLI_none) +  [-errX, errX],...
    [mean(meanNLI_nat), mean(meanNLI_nat)],'errX',fig4,'g','-','none')

addLineToAxis([0 1],[0 1],'unity',fig4,'k','--','none')
[p, ~] = signrank(meanNLI_none,meanNLI_nat)
figID = ['MixSurMod_no_nat',figTag];
makeAxisStruct(fig4,figID ,'RFSurroundFigs')


figure; clf;
fig5=gca;
set(fig5,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig5,'XLabel'),'String','NLI no surround')
set(get(fig5,'YLabel'),'String','NLI random surround')
addLineToAxis(meanNLI_none,meanNLI_mix,'data',fig5,'r','none','o')

errY = std(meanNLI_mix) ./ sqrt(length(meanNLI_mix));
errX = std(meanNLI_none) ./ sqrt(length(meanNLI_none));

addLineToAxis(mean(meanNLI_none),mean(meanNLI_mix),'mean',fig5,'r','none','s')

addLineToAxis([mean(meanNLI_none), mean(meanNLI_none)],...
    mean(meanNLI_mix) + [-errY, errY],'errY',fig5,'r','-','none')
addLineToAxis(mean(meanNLI_none) +  [-errX, errX],...
    [mean(meanNLI_mix), mean(meanNLI_mix)],'errX',fig5,'r','-','none')

addLineToAxis([0 1],[0 1],'unity',fig5,'k','--','none')
[p, ~] = signrank(meanNLI_none,meanNLI_mix)
figID = ['MixSurMod_no_mix',figTag];
makeAxisStruct(fig5,figID ,'RFSurroundFigs')

figure; clf;
fig6=gca;
set(fig6,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig6,'XLabel'),'String','NLI natural surround')
set(get(fig6,'YLabel'),'String','NLI random surround')
addLineToAxis(meanNLI_nat,meanNLI_mix,'data',fig6,'k','none','o')

errY = std(meanNLI_mix) ./ sqrt(length(meanNLI_mix));
errX = std(meanNLI_nat) ./ sqrt(length(meanNLI_nat));

addLineToAxis(mean(meanNLI_nat),mean(meanNLI_mix),'mean',fig6,'k','none','s')

addLineToAxis([mean(meanNLI_nat), mean(meanNLI_nat)],...
    mean(meanNLI_mix) + [-errY, errY],'errY',fig6,'k','-','none')
addLineToAxis(mean(meanNLI_nat) +  [-errX, errX],...
    [mean(meanNLI_mix), mean(meanNLI_mix)],'errX',fig6,'k','-','none')

addLineToAxis([0 1],[0 1],'unity',fig6,'k','--','none')
[p, ~] = signrank(meanNLI_nat,meanNLI_mix)
figID = ['MixSurMod_nat_mix',figTag];
makeAxisStruct(fig6,figID ,'RFSurroundFigs')


%%
noBins = 40;
binEdges = linspace(-1,1,noBins+1);
binCtrs = binEdges(1:end-1) + diff(binEdges);

results.noDiff = [];
results.noNLI= [];

results.natDiff = [];
results.natNLI= [];

results.mixDiff = [];
results.mixNLI = [];

% figure(1); clf;
% figure(2); clf;
ct = 0;
for currentImageIndex = 1:30
    if currentImageIndex == 12
        continue
    end
    ct = ct + 1;
    
    load(['ImageDiscModel_',num2str(currentImageIndex),'_20170718.mat'])

    intDiff_no = (stats.centerMean) ./ stats.imageMean;
    intDiff_nat = (stats.centerMean - stats.surroundMean) ./ stats.imageMean;
    intDiff_mix = (stats.centerMean - stats.surroundMean_Mix) ./ stats.imageMean;

    %calc NLIs:
    NLI_noSurround = (response.Image - response.Disc) ./...
        (response.Image + response.Disc);
    
    NLI_naturalSurround = (response.ImageMatchedSurround - response.DiscMatchedSurround) ./...
        (response.ImageMatchedSurround + response.DiscMatchedSurround);

    NLI_mixedSurround = (response.ImageMixedSurround - response.DiscMixedSurround) ./...
        (response.ImageMixedSurround + response.DiscMixedSurround);
% 
% NLI_naturalSurround = (response.ImageMatchedSurround - response.DiscMatchedSurround);
% 
%     NLI_mixedSurround = (response.ImageMixedSurround - response.DiscMixedSurround);
%     
    results.noDiff = cat(1,results.noDiff,intDiff_no);
    results.noNLI = cat(1,results.noNLI,NLI_noSurround);
    
    results.natDiff = cat(1,results.natDiff,intDiff_nat);
    results.natNLI = cat(1,results.natNLI,NLI_naturalSurround);
    
    results.mixDiff = cat(1,results.mixDiff,intDiff_mix);
    results.mixNLI = cat(1,results.mixNLI,NLI_mixedSurround);
    
end


figure; clf;
fig6=gca;
set(fig6,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig6,'XLabel'),'String','Ic - Is')
set(get(fig6,'YLabel'),'String','NLI')

results.natDiff = results.natDiff(:);
results.natNLI = results.natNLI(:);
keepInds = find(~isnan(results.natNLI));
results.natDiff = results.natDiff(keepInds);
results.natNLI = results.natNLI(keepInds);
binAndPlotEquallyPopulatedBins(results.natDiff(:),results.natNLI(:),40, fig6 ,'g');

results.mixDiff = results.mixDiff(:);
results.mixNLI = results.mixNLI(:);
keepInds = find(~isnan(results.mixNLI));
results.mixDiff = results.mixDiff(keepInds);
results.mixNLI = results.mixNLI(keepInds);
binAndPlotEquallyPopulatedBins(results.mixDiff(:),results.mixNLI(:),40, fig6, 'r');

figID = ['MixSur_NLIvsDiff'];
% makeAxisStruct(fig6,figID ,'RFSurroundFigs')

figure; clf;
fig7=gca;
set(fig7,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig7,'XLabel'),'String','Ic - Is')
set(get(fig7,'YLabel'),'String','NLI')

results.noDiff = results.noDiff(:);
results.noNLI = results.noNLI(:);
keepInds = find(~isnan(results.noNLI));
results.noDiff = results.noDiff(keepInds);
results.noNLI = results.noNLI(keepInds);
binAndPlotEquallyPopulatedBins(results.noDiff(:),results.noNLI(:),40, fig7 ,'k');


figID = ['MixSur_NLIvsDiff'];

