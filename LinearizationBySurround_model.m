clear all; clc;
load('NaturalImageFlashLibrary_101716.mat')
imageNames = fieldnames(imageData);
resourcesDir = '~/Documents/MATLAB/turner-package/resources';
patchSize = 120; %pixels
MicronsPerPixel = 6.6;
patchesPerImage = 1e3;
imageIndex = 1; %up to 20
cellPolarity = -1;
rng(1); %set seed

RFmodel = ThreeLayerReceptiveFieldModel;
RFmodel.MicronsPerPixel = MicronsPerPixel;
RFmodel.makeRfComponents(patchSize);

ImageX = 1536; ImageY = 1024;

respFieldNames = {'CenterOnly_LN','CenterOnly_NonlinearSubunits',...
    'CenterSurround_LN',...
    'CenterSurround_NonlinearCenterPlusNonlinearSurround',...
    'CenterSurround_SharedNonlinearity',...
    'CenterSurround_LNLN',...
    'CenterSurround_NonlinearCenterPlusIndependentLinearSurround'};
ImageID = imageNames{imageIndex};
fileId=fopen([resourcesDir, '/VHsubsample_20160105', '/', ImageID,'.iml'],'rb','ieee-be');
img = fread(fileId, [1536,1024], 'uint16');
contrastImg = img - mean(img(:));
contrastImg = contrastImg ./ max(contrastImg(:));

location = zeros(patchesPerImage,2);
tempResponse = struct;
for cc = 1:length(respFieldNames)
    tempResponse.(respFieldNames{cc}) = zeros(patchesPerImage,1);
end
for rr = 1:patchesPerImage
     % choose location randomly
    x = round(patchSize/2 + (ImageX - patchSize)*rand);
    y = round(patchSize/2 + (ImageY - patchSize)*rand);
    location(rr,:) = [x, y];

    % get patch
    ContrastPatch = cellPolarity .* contrastImg(x-patchSize/2+1:x+patchSize/2,y-patchSize/2+1:y+patchSize/2);

    % get RF model responses
    responseStructure = RFmodel.getResponse(ContrastPatch);
    for cc = 1:length(respFieldNames)
        tempResponse.(respFieldNames{cc})(rr) = responseStructure.(respFieldNames{cc});
    end
end

% pull relatively nonlinear patches with similar response strengths...
targetLinearResponse = 0.015; targetSubunitResponse = 0.04;
% % tempDistance = sqrt((tempResponse.CenterOnly_LN - targetLinearResponse).^2 + ...
% %     (tempResponse.CenterOnly_NonlinearSubunits - targetSubunitResponse).^2);

tempDistance = sqrt((tempResponse.CenterOnly_NonlinearSubunits - targetSubunitResponse).^2);

[val indsTemp] = sort(abs(tempDistance));
pullInds = indsTemp(1:20);

figure(1); clf;
subplot(131); hold on;
limUp = max([tempResponse.CenterOnly_NonlinearSubunits,tempResponse.CenterOnly_LN]);
plot([0 limUp],[0 limUp],'k--')
plot(tempResponse.CenterOnly_NonlinearSubunits,tempResponse.CenterOnly_LN,...
    'Marker','o','LineStyle','none','Color','k');

plot(tempResponse.CenterOnly_NonlinearSubunits(pullInds),tempResponse.CenterOnly_LN(pullInds),...
    'Marker','.','LineStyle','none','Color','r');
title('Center only'); xlabel('Nonlinear subunits'); ylabel('LN')

subplot(132); hold on;
limUp = max([tempResponse.CenterSurround_SharedNonlinearity,tempResponse.CenterSurround_LN]);
plot([0 limUp],[0 limUp],'k--')
plot(tempResponse.CenterSurround_SharedNonlinearity,tempResponse.CenterSurround_LN,...
    'Marker','o','LineStyle','none','Color','k');
title('Center + surround'); xlabel('Subunits with surround'); ylabel('LN C+S')

tempCenter = tempResponse.CenterOnly_NonlinearSubunits(pullInds);
tempShared = tempResponse.CenterSurround_SharedNonlinearity(pullInds);
tempIndep = tempResponse.CenterSurround_NonlinearCenterPlusIndependentLinearSurround(pullInds);

%sort by shared NL model:
[val, ind] = sort(tempShared,'descend');

figure(2); clf; hold on;
plot(1:20,tempCenter(ind),'ko')
plot(1:20,tempShared(ind),'ro')
plot(1:20,tempIndep(ind),'bo')

%% get parameters of RF models s.t. surround strength is == across models

corrFieldNames = {'imageContrast','LNcenter','subunitCenter','LNCenterSurround',...
    'SubunitCenterPlusLinearSurround','SubunitCenterPlusSubunitSurround','LinearCenterPlusNonlinearSurround'};

patchSize = 120; %pixels
MicronsPerPixel = 6.6;
RFmodel = ThreeLayerReceptiveFieldModel;
RFmodel.MicronsPerPixel = MicronsPerPixel;
RFmodel.makeRfComponents(patchSize);

% expanding spots...
spotSize = [5:5:patchSize];
[rr, cc] = meshgrid(1:patchSize,1:patchSize);
response.LNCenterSurround = [];
response.SubunitCenterPlusLinearSurround = [];
response.SubunitCenterPlusSubunitSurround = [];
response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround = [];
for ss = 1:length(spotSize) %get responses to each spot
    currentRadius = spotSize(ss)/2;
    spotBinary = double(sqrt((rr-(patchSize/2)).^2+(cc-(patchSize/2)).^2)<=currentRadius);
    responseStructure = RFmodel.getResponse(spotBinary);
    
    response.LNCenterSurround(ss) = responseStructure.CenterSurround_LN;
    response.SubunitCenterPlusLinearSurround(ss) = responseStructure.CenterSurround_SharedNonlinearity;
    response.SubunitCenterPlusSubunitSurround(ss) = responseStructure.CenterSurround_LNLN;
    response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround(ss) = responseStructure.CenterSurround_NonlinearCenterPlusIndependentLinearSurround;
end
spotSize = spotSize .* MicronsPerPixel;
figure(2); clf; hold on
plot(spotSize,response.LNCenterSurround ./ max(response.LNCenterSurround),'b-')
plot(spotSize,response.SubunitCenterPlusLinearSurround ./ max(response.SubunitCenterPlusLinearSurround),'r-o')
plot(spotSize,response.SubunitCenterPlusSubunitSurround ./ max(response.SubunitCenterPlusSubunitSurround),'k-')
plot(spotSize,response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround ./ ...
    max(response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround),'g-')
legend('LN','SharedNL','LNLN','IndepNL')
%%
% gratings...
barSizes = [5 10 20 30 40 50 60 70 80 90 100,...
    120 140 160 180 200 220 240 260 280 300 320];
response.LNCenterSurround = [];
response.SubunitCenterPlusLinearSurround = [];
response.SubunitCenterPlusSubunitSurround = [];
response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround = [];
for bb = 1:length(barSizes);
    barSize = barSizes(bb); %microns
    sf = 1 / (2*barSize/MicronsPerPixel);
    tempXX = sin(2*pi.*sf*(-patchSize/2 + 1: patchSize/2));
    newStim = repmat(tempXX,patchSize,1);
    responseStructure = RFmodel.getResponse(newStim);
    
    response.LNCenterSurround(bb) = responseStructure.CenterSurround_LN;
    response.SubunitCenterPlusLinearSurround(bb) = responseStructure.CenterSurround_SharedNonlinearity;
    response.SubunitCenterPlusSubunitSurround(bb) = responseStructure.CenterSurround_LNLN;
    response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround(bb) = responseStructure.CenterSurround_NonlinearCenterPlusIndependentLinearSurround;
end

figure(2); clf; hold on
plot(barSizes,response.LNCenterSurround ./ max(response.LNCenterSurround),'b-')
plot(barSizes,response.SubunitCenterPlusLinearSurround ./ max(response.SubunitCenterPlusLinearSurround),'r-o')
plot(barSizes,response.SubunitCenterPlusSubunitSurround ./ max(response.SubunitCenterPlusSubunitSurround),'k-')
plot(barSizes,response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround ./ ...
    max(response.CenterSurround_NonlinearCenterPlusIndependentLinearSurround),'g-')

legend('LN','SharedNL','LNLN','IndepNL')
