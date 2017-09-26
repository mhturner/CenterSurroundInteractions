%% RF components
clear all; close all; clc;

%sub sur: 0.10, 0.72, 1.1
%ind sur: 0.10, 0.70, 1.04
SubunitSurroundWeight = 0.72;     % relative to center integral
contrastPolarity = -1;

%filter sizes
filterSize = 300;                
subunitSigma = 3;
subunitSurroundSigma = 75;
centerSigma = 20;

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
    end
end

% normalize components
SubunitFilter = SubunitFilter / sum(SubunitFilter(:));
SubunitSurroundFilter = SubunitSurroundFilter / sum(SubunitSurroundFilter(:));
RFCenter = RFCenter / sum(RFCenter(:));

SubunitWithSurroundFilter = SubunitFilter - SubunitSurroundWeight .* SubunitSurroundFilter;

%get weighting of each subunit output by center
subunitWeightings = RFCenter(SubunitIndices);
subunitWeightings = subunitWeightings ./ sum(subunitWeightings);

%% Make grating stimulus
f = 12;
dim = filterSize;

per = dim / f;
disp(per);

tt = linspace(0,dim,dim);
a = sin(2 * pi .* tt .* f);
a(a>=0) = 0.25;
a(a<0) = -0.25
imageMat = repmat(a,dim,1);

figure(1); clf; 
subplot(3,1,1);
plot(tt,a,'k')
subplot(3,1,[2,3])
imagesc(imageMat); colormap(gray);


%% responses of RF models
tic;

response.Center_LN = [];
response.Center_subunit = [];

response.Center_SharedSurround = [];

bg = linspace(-0.5,0.5,21);

for pp = 1:length(bg)
        
    % Grating responses:
    newImagePatch_subunit = contrastPolarity.*conv2(imageMat + bg(pp),SubunitFilter,'same');
    newImagePatch_subSurround = contrastPolarity.*conv2(imageMat + bg(pp),SubunitWithSurroundFilter,'same');
    
    % activation of each subunit
    subunitActivations = newImagePatch_subunit(SubunitIndices);
    subunitOutputs = subunitActivations;
    
    %Center 1 - linear subunit center and output NL:
    response.Center_LN(pp) = max(sum(subunitOutputs .* subunitWeightings),0); %post-summation rectification
    %Center 2 - nonlinear subunit center:
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_subunit(pp) = sum(subunitOutputs .* subunitWeightings);
 
    %CenterSurround 4 - Shared NL, i.e. subunits have surrounds:
    % activation of each subunit
    subunitActivations = newImagePatch_subSurround(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_SharedSurround(pp) = sum(subunitOutputs .* subunitWeightings);
    
    
    % Mean responses:
    newImagePatch_subunit = contrastPolarity.*conv2(zeros(size(imageMat)) + bg(pp),SubunitFilter,'same');
    newImagePatch_subSurround = contrastPolarity.*conv2(zeros(size(imageMat)) + bg(pp),SubunitWithSurroundFilter,'same');
    
    % activation of each subunit
    subunitActivations = newImagePatch_subunit(SubunitIndices);
    subunitOutputs = subunitActivations;
    
    %Center 1 - linear subunit center and output NL:
    response.Center_LN_mean(pp) = max(sum(subunitOutputs .* subunitWeightings),0); %post-summation rectification
    %Center 2 - nonlinear subunit center:
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_subunit_mean(pp) = sum(subunitOutputs .* subunitWeightings);

    %CenterSurround 4 - Shared NL, i.e. subunits have surrounds:
    % activation of each subunit
    subunitActivations = newImagePatch_subSurround(SubunitIndices);
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0; %threshold each subunit
    response.Center_SharedSurround_mean(pp) = sum(subunitOutputs .* subunitWeightings);
   
end
toc;
%%
figure; clf;
fig3=gca;
set(fig3,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig3,'XLabel'),'String','Mean Intensity')
set(get(fig3,'YLabel'),'String','Response (norm.)')

addLineToAxis(bg,response.Center_subunit./max(response.Center_subunit),...
    'centerGrating',fig3,'k','-','none')
addLineToAxis(bg,response.Center_subunit_mean./max(response.Center_subunit_mean),...
    'centerMean',fig3,'k','--','none')

addLineToAxis(bg,response.Center_SharedSurround./max(response.Center_SharedSurround),...
    'csGrating',fig3,'g','-','none')
addLineToAxis(bg,response.Center_SharedSurround_mean./max(response.Center_SharedSurround_mean),...
    'csMean',fig3,'g','--','none')

makeAxisStruct(fig3,'GratingModelResps' ,'RFSurroundFigs')


figure; clf;
fig4=gca;
set(fig4,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig4,'XLabel'),'String','Mean Intensity')
set(get(fig4,'YLabel'),'String','Response (norm.)')
addLineToAxis(tt,a + bg(1),'stim',fig4,'k','-','none')
makeAxisStruct(fig4,'GratingModelStims1' ,'RFSurroundFigs')


figure; clf;
fig5=gca;
set(fig5,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig5,'XLabel'),'String','Mean Intensity')
set(get(fig5,'YLabel'),'String','Response (norm.)')
addLineToAxis(tt,a + bg(11),'stim',fig5,'k','-','none')
makeAxisStruct(fig5,'GratingModelStims2' ,'RFSurroundFigs')

figure; clf;
fig6=gca;
set(fig6,'XScale','linear','YScale','linear')
set(0, 'DefaultAxesFontSize', 12)
set(get(fig6,'XLabel'),'String','Mean Intensity')
set(get(fig6,'YLabel'),'String','Response (norm.)')
addLineToAxis(tt,a + bg(end),'stim',fig6,'k','-','none')
makeAxisStruct(fig6,'GratingModelStims3' ,'RFSurroundFigs')


