clear all; close all; clc;
tic;
%eg stims: 1, 50, 35, 10
stimInd = 35; %stim to use
micronsPerPixel = 6.6; %6.6

contrastPolarity = -1;

%RF components:
%sub sur: 0.10, 0.72, 1.1
%ind sur: 0.10, 0.70, 1.04
SubunitSurroundWeight = 0.72;     % relative to center integral
SurroundWeight = 0.70;           % relative to center integral

filterSize_um = 700;                % size of patch (um)
subunitSigma_um = 10;
subunitSurroundSigma_um = 150;
centerSigma_um = 40;
surroundSigma_um = 150;

figDir = '~/Dropbox/RiekeLab/Analysis/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs


resources_dir = '~/Dropbox/RiekeLab/Analysis/MATLAB/turner-package/resources/';
IMAGES_DIR_VH = '~/Dropbox/RiekeLab/Analysis/MATLAB/MHT-analysis/resources/vanhateren_iml/';
load([resources_dir, 'dovesFEMstims_20160422.mat'])

%image:
fileName = [FEMdata(stimInd).ImageName(1:8), '_',num2str(FEMdata(stimInd).SubjectIndex),'.mat'];
f1=fopen([IMAGES_DIR_VH, FEMdata(stimInd).ImageName],'rb','ieee-be');
w=1536;h=1024;
my_image=fread(f1,[w,h],'uint16');
my_image = my_image';

%trim it square
my_image = my_image(1:1024,1:1024);
w = 1024; h = 1024;


my_image_nomean = (my_image - mean(my_image(:))) ./ mean(my_image(:));

my_image_nomean  = contrastPolarity .* my_image_nomean;

%convert to pixels:
filterSize = round(filterSize_um / micronsPerPixel);                % size of patch (um)
subunitSigma = round(subunitSigma_um / micronsPerPixel);
subunitSurroundSigma = round(subunitSurroundSigma_um / micronsPerPixel);
centerSigma = round(centerSigma_um / micronsPerPixel);
surroundSigma = round(surroundSigma_um / micronsPerPixel);


TempFilter = zeros(filterSize, filterSize);
SubunitLocations = find(rem([1:filterSize], 2*subunitSigma) == 0);
for x = 1:length(SubunitLocations)
    TempFilter(SubunitLocations(x), SubunitLocations) = 1;
end
SubunitIndices = find(TempFilter > 0);

% Filters for subunit, subunit with surround, center, surround
for x = 1:filterSize
    for y = 1:filterSize
        SubunitFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (subunitSigma^2)));
        SubunitSurroundFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (subunitSurroundSigma^2)));
        
        RFCenter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (centerSigma^2)));
        SurroundFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (surroundSigma^2))); 
    end
end

% normalize components
SubunitFilter = SubunitFilter / sum(SubunitFilter(:));
SubunitSurroundFilter = SubunitSurroundFilter / sum(SubunitSurroundFilter(:));
SurroundFilter = SurroundFilter  / sum(SurroundFilter(:));
RFCenter = RFCenter / sum(RFCenter(:));

SubunitWithSurroundFilter = SubunitFilter - SubunitSurroundWeight .* SubunitSurroundFilter;
SurroundFilter = SurroundFilter .* SurroundWeight;

%get weighting of each subunit output by center
subunitWeightings = RFCenter(SubunitIndices);
subunitWeightings = subunitWeightings ./ sum(subunitWeightings);

%get filter-convolved image(s)
convolved_surround = conv2(my_image_nomean,SurroundFilter,'same');

convolved_subunit = conv2(my_image_nomean,SubunitFilter,'same');
convolved_subSurround = conv2(my_image_nomean,SubunitWithSurroundFilter,'same');
toc;


figure(1); clf; imagesc(my_image); colormap(gray)
%% responses to image patches
tic;

response.Center_LN = [];
response.Center_subunit = [];

response.Center_PostNLSurround = [];
response.Center_SharedSurround = [];

xx = (filterSize/2) : (w-filterSize/2);
yy = (filterSize/2) : (h-filterSize/2);

[X, Y] = meshgrid(xx,yy);
patchLocation = [];
patchLocation(:,1) = X(:);
patchLocation(:,2) = Y(:);

noPatches = length(patchLocation);

for pp = 1:noPatches
    x = patchLocation(pp,1);
    y = patchLocation(pp,2);
    
    newImagePatch_surround = convolved_surround(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    newImagePatch_subunit = convolved_subunit(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    newImagePatch_subSurround = convolved_subSurround(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    
    
    % activation of each subunit
    subunitActivations = newImagePatch_subunit(SubunitIndices);
    subunitOutputs = subunitActivations;
    
    %Center 1 - linear subunit center and output NL:
    response.Center_LN(y,x) = max(sum(subunitOutputs .* subunitWeightings),0); %post-summation rectification
    %Center 2 - nonlinear subunit center:
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_subunit(y,x) = sum(subunitOutputs .* subunitWeightings);
    
    %CenterSurround 3 - nonlinear subunit center and post-NL surround:
    temp = sum(subunitOutputs .* subunitWeightings) - (newImagePatch_surround(filterSize/2,filterSize/2));
    response.Center_PostNLSurround(y,x) = max(temp,0); %output nonlinearity
    
    %CenterSurround 3.5 - linear subunit center and surround. Output NL (~DoG):
    temp = sum(subunitActivations .* subunitWeightings) - (newImagePatch_surround(filterSize/2,filterSize/2));
    response.DoG(y,x) = max(temp,0); %output nonlinearity
    
    %CenterSurround 4 - Shared NL, i.e. subunits have surrounds:
    % activation of each subunit
    subunitActivations = newImagePatch_subSurround(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_SharedSurround(y,x) = sum(subunitOutputs .* subunitWeightings);
   
end
toc;

%%
buf = filterSize/2;

fh = figure(9); clf;
imagesc(my_image(buf+1:end-buf,buf+1:end-buf)); colormap(gray);
brighten(0.6)
xlim([0,500]); ylim([0 500]);
axis image; axis off;
drawnow;
figID = ['patch_fig_',num2str(stimInd)];
print(fh,[figDir,figID],'-depsc')

fh = figure(8);
imagesc(response.Center_subunit(buf+1:end,buf+1:end)); colormap(gray);
caxis([0,1.1*max(response.Center_subunit(:))]); colorbar;
brighten(0.6)
xlim([0,500]); ylim([0 500]);
axis image; axis off;
drawnow;
figID = ['patch_cbar_',num2str(stimInd)];
print(fh,[figDir,figID],'-depsc')


fh = figure(10);
imagesc(response.Center_subunit(buf+1:end,buf+1:end)); colormap(gray);
caxis([0,1.1*max(response.Center_subunit(:))]);
brighten(0.6)
xlim([0,500]); ylim([0 500]);
axis image; axis off;
drawnow;
figID = ['patch_sub_',num2str(stimInd)];
print(fh,[figDir,figID],'-depsc')

% fh = figure(11);
% imagesc(response.Center_PostNLSurround(buf+1:end,buf+1:end)); colormap(gray);
% caxis([0,1.1*max(response.Center_subunit(:))]);
% brighten(0.6)
% axis image; axis off;
% drawnow;
% figID = ['patch_postNL_',num2str(stimInd)];
% print(fh,[figDir,figID],'-depsc')

fh = figure(12);
imagesc(response.Center_SharedSurround(buf+1:end,buf+1:end)); colormap(gray);
caxis([0,1.1*max(response.Center_subunit(:))]);
brighten(0.6)
xlim([0,500]); ylim([0 500]);
axis image; axis off;
drawnow;
figID = ['patch_sharedNL_',num2str(stimInd)];
print(fh,[figDir,figID],'-depsc')

fh = figure(13);
imagesc(response.Center_LN(buf+1:end,buf+1:end)); colormap(gray);
caxis([0,1.1*max(response.Center_subunit(:))]);
brighten(0.6)
axis image; axis off;
drawnow;
figID = ['patch_centerLN_',num2str(stimInd)];
print(fh,[figDir,figID],'-depsc')

fh = figure(14);
imagesc(response.DoG(buf+1:end,buf+1:end)); colormap(gray);
caxis([0,1.1*max(response.Center_subunit(:))]);
brighten(0.6)
axis image; axis off;
drawnow;
figID = ['patch_DoG_',num2str(stimInd)];
print(fh,[figDir,figID],'-depsc')

%% power spectra:
trimSize = 970;

figure; clf; fig1=gca;
set(fig1,'XScale','log','YScale','log')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig1,'XLabel'),'String','Spatial frequency')
set(get(fig1,'YLabel'),'String','Power')

tempP = getPowerSpectrum(my_image(1:trimSize,1:trimSize),6.6);
addLineToAxis(tempP.f,tempP.p./tempP.p(2),'Image',fig1,'k','-','none')

tempP = getPowerSpectrum(response.Center_subunit(1:trimSize,1:trimSize),6.6);
addLineToAxis(tempP.f,tempP.p./tempP.p(2),'CenterSubunit',fig1,'b','-','none')

tempP = getPowerSpectrum(response.DoG(1:trimSize,1:trimSize),6.6);
addLineToAxis(tempP.f,tempP.p./tempP.p(2),'DoG',fig1,'g','-','none')

tempP = getPowerSpectrum(response.Center_SharedSurround(1:trimSize,1:trimSize),6.6);
addLineToAxis(tempP.f,tempP.p./tempP.p(2),'SharedSurround',fig1,'r','-','none')
%% Exp spots...
figure(2); clf;  set(gcf, 'WindowStyle', 'docked')
subplot(221); imagesc(RFCenter); colormap(gray); colorbar; title('Center')
subplot(222); imagesc(SurroundFilter); colormap(gray); colorbar; title('Surround')
subplot(223); imagesc(SubunitFilter); colormap(gray); colorbar; title('Subunit')
subplot(224); imagesc(SubunitWithSurroundFilter); colormap(gray); colorbar; title('Subunit w/ surround')

% response to spots
spotSize = 0:10:600; %diameter of spot

Response_model1 = [];
Response_model2 = [];

[rr, cc] = meshgrid(1:filterSize,1:filterSize);
for ss = 1:length(spotSize) %get responses to each spot
    currentRadius = (spotSize(ss)/2)/ micronsPerPixel; %convert to pixel
    spotBinary = double(sqrt((rr-(filterSize/2)).^2+(cc-(filterSize/2)).^2)<=currentRadius);
    
    %Model 1: Subunits combine, then linear surround:
    convolved_Subunit = conv2(spotBinary, SubunitFilter, 'same');  
    convolved_Surround = conv2(spotBinary, SurroundFilter, 'same');
        % activation of each subunit
    subunitActivations = convolved_Subunit(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
        %Add the linear surround to the subunit center output:
    temp = sum(subunitOutputs .* subunitWeightings) - max(convolved_Surround(filterSize/2,filterSize/2),0);
    Response_model1(ss) = max(temp,0); %output nonlinearity
    
    %Model 2: Subunits have surrounds:
    convolved_SubunitWithSurround = conv2(spotBinary, SubunitWithSurroundFilter, 'same');
        % activation of each subunit
    subunitActivations = convolved_SubunitWithSurround(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    Response_model2(ss) = sum(subunitOutputs .* subunitWeightings);
end

figure(4);  clf; set(gcf, 'WindowStyle', 'docked')
plot(spotSize,Response_model1,'k-'); hold on;
plot(spotSize,Response_model2,'b-')
legend('Post NL surround','Shared NL')

