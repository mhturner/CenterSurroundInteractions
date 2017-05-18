clear all; close all; clc;

stimInd = 1; %image to use
noPatches = 5000;
micronsPerPixel = 6.6;

%RF components:
SubunitSurroundWeight = 0.72;     % relative to center integral
SurroundWeight = 0.7;           % relative to center integral

filterSize_um = 700;                % size of patch (um)
subunitSigma_um = 10;
subunitSurroundSigma_um = 150;
centerSigma_um = 40;
surroundSigma_um = 150;


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
    Response_model1(ss) = sum(subunitOutputs .* subunitWeightings) - convolved_Surround(filterSize/2,filterSize/2);

    
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

%% responses to image patches
tic;

response.Center_LN = [];
response.Center_subunit = [];

response.Center_PostNLSurround = [];
response.Center_SharedSurround = [];

patchLocation = [];
for pp = 1:noPatches
    x = round(filterSize/2 + (w - filterSize)*rand);
    y = round(filterSize/2 + (h - filterSize)*rand);
    newImagePatch = my_image_nomean(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
    patchLocation(pp,:) = [x, y];
    
    % non-surround subunits and linear surround
    convolved_Subunit = conv2(newImagePatch, SubunitFilter, 'same');  
    convolved_Surround = conv2(newImagePatch, SurroundFilter, 'same');
    % activation of each subunit
    subunitActivations = convolved_Subunit(SubunitIndices);
    subunitOutputs = subunitActivations;
    
    %Center 1 - linear subunit center and output NL:
    response.Center_LN(pp) = max(sum(subunitOutputs .* subunitWeightings),0); %post-summation rectification
    %Center 2 - nonlinear subunit center:
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_subunit(pp) = sum(subunitOutputs .* subunitWeightings);
    
    %CenterSurround 3 - nonlinear subunit center and post-NL surround:
    response.Center_PostNLSurround(pp) = sum(subunitOutputs .* subunitWeightings) - max(convolved_Surround(filterSize/2,filterSize/2),0);
    
    %CenterSurround 4 - Shared NL, i.e. subunits have surrounds:
    convolved_SubunitWithSurround = conv2(newImagePatch, SubunitWithSurroundFilter, 'same');
    % activation of each subunit
    subunitActivations = convolved_SubunitWithSurround(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_SharedSurround(pp) = sum(subunitOutputs .* subunitWeightings);
   
end
toc;

%%
figure(5); clf;  set(gcf, 'WindowStyle', 'docked')
subplot(221);
plot(response.Center_subunit,response.Center_LN,'ko'); hold on;
plot([0 1],[0 1],'k--');
xlabel('Subunit center'); ylabel('LN center');

subplot(222);
plot(response.Center_SharedSurround, response.Center_PostNLSurround, 'ko'); hold on
plot([0 0.8],[0 0.8],'k--')
xlabel('Shared NL'); ylabel('Post-NL surround');

% diff between surround models:
diff = response.Center_SharedSurround - response.Center_PostNLSurround;
[~, ind] = sort(diff);

cMat = jet(noPatches);

figure(3); clf;
imagesc(my_image); axis image; axis equal; axis off; colormap(gray);
brighten(0.6) %brighten for display purposes
hold on;
scatter(patchLocation(ind,1),patchLocation(ind,2),40*(1:noPatches)./(noPatches),cMat);


%%

diff = response.Center_SharedSurround - response.Center_PostNLSurround;
[val, ind] = sort(diff,'descend');
lookInds = ind(1:10);
lowInds = ind(end-10:end);

   figure(10); clf;
for pp = 1:10
   x = patchLocation(lookInds(pp),1);
   y = patchLocation(lookInds(pp),2);
   newImagePatch = my_image(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
   subplot(4,5,pp);
   imagesc(newImagePatch); colormap(gray); axis image; axis equal; axis off;
end

for pp = 1:10
   x = patchLocation(ind(end-pp),1);
   y = patchLocation(ind(end-pp),2);
   newImagePatch = my_image(y - filterSize/2 + 1 : y + filterSize/2, x - filterSize/2 + 1 : x + filterSize/2);
   subplot(4,5,10+pp);
   imagesc(newImagePatch); colormap(gray); axis image; axis equal; axis off;
end


