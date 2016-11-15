clear all; clc;
load('NaturalImageFlashLibrary_101716.mat')
imageNames = fieldnames(imageData);
resourcesDir = '~/Documents/MATLAB/turner-package/resources';
patchSize = 100; %pixels
MicronsPerPixel = 6.6;
patchesPerImage = 1e2;
noImages = 5; %up to 20

RFmodel = ThreeLayerReceptiveFieldModel;
RFmodel.MicronsPerPixel = MicronsPerPixel;
RFmodel.makeRfComponents(patchSize);

ImageX = 1536; ImageY = 1024;

colors = pmkmp(20);
respFieldNames = {'LNcenter','subunitCenter','LNCenterSurround',...
    'SubunitCenterPlusLinearSurround','SubunitCenterPlusSubunitSurround','LinearCenterPlusNonlinearSurround',...
    'DivLL','DivLN','DivNL','DivNN'};
figure(1); clf;
for ii = 3%1:noImages
    disp(num2str(ii))
    ImageID = imageNames{ii};
    fileId=fopen([resourcesDir, '/VHsubsample_20160105', '/', ImageID,'.iml'],'rb','ieee-be');
    img = fread(fileId, [1536,1024], 'uint16');
    contrastImg = img - mean(img(:));
    contrastImg = contrastImg ./ max(contrastImg(:));

    noBins = 2.*patchesPerImage; %from no. image patches to show
    [N, edges, bin] = histcounts(imageData.(ImageID).responseDifferences,noBins);
    populatedBins = unique(bin);
    %pluck one patch from each bin
    pullInds = arrayfun(@(b) find(b == bin,1),populatedBins);
    location = imageData.(ImageID).location(pullInds,:);
    
    
    tempResponse = struct;
    for cc = 1:length(respFieldNames)
        tempResponse.(respFieldNames{cc}) = zeros(patchesPerImage,1);
    end
    for rr = 1:patchesPerImage
         % choose location
%         x = round(patchSize/2 + (ImageX - patchSize)*rand);
%         y = round(patchSize/2 + (ImageY - patchSize)*rand);
%         location(rr,:) = [x, y];

        x = location(rr,1); y = location(rr,2);

        % get patch
        ContrastPatch = -contrastImg(x-patchSize/2+1:x+patchSize/2,y-patchSize/2+1:y+patchSize/2);
        
        % get RF model responses
        responseStructure = RFmodel.getResponse(ContrastPatch);
        tempResponse.LNcenter(rr) = responseStructure.CenterOnly.LN;
        tempResponse.subunitCenter(rr) = responseStructure.CenterOnly.NonlinearSubunits;
        
        tempResponse.LNCenterSurround(rr) = responseStructure.CenterSurround.LN;
        tempResponse.SubunitCenterPlusLinearSurround(rr) = responseStructure.CenterSurround.SharedNonlinearity;
        tempResponse.SubunitCenterPlusSubunitSurround(rr) = responseStructure.CenterSurround.LNLNSamePolarity;
        tempResponse.LinearCenterPlusNonlinearSurround(rr) = responseStructure.CenterSurround.LinearCenterPlusNonlinearSurround;
        
        tempResponse.DivLL(rr) = responseStructure.DivisiveSurround.linC_linS;
        tempResponse.DivLN(rr) = responseStructure.DivisiveSurround.linC_nlS;
        tempResponse.DivNL(rr) = responseStructure.DivisiveSurround.nlC_linS;
        tempResponse.DivNN(rr) = responseStructure.DivisiveSurround.nlC_nlS;
    end
    figure(1);
    subplot(131); hold on;
    limUp = max([tempResponse.subunitCenter,tempResponse.LNcenter]);
    plot([0 limUp],[0 limUp],'k--')
    plot(tempResponse.subunitCenter,tempResponse.LNcenter,...
        'Marker','o','LineStyle','none','Color',colors(ii,:));
    
    subplot(132); hold on;
    limUp = max([tempResponse.SubunitCenterPlusLinearSurround,tempResponse.LNCenterSurround]);
    plot([0 limUp],[0 limUp],'k--')
    plot(tempResponse.SubunitCenterPlusLinearSurround,tempResponse.LNcenter,...
        'Marker','o','LineStyle','none','Color',colors(ii,:));
    
    subplot(133); hold on;
    limUp = max([tempResponse.SubunitCenterPlusSubunitSurround,tempResponse.LinearCenterPlusNonlinearSurround]);
    plot([0 limUp],[0 limUp],'k--')
    plot(tempResponse.SubunitCenterPlusSubunitSurround,tempResponse.LNcenter,...
        'Marker','o','LineStyle','none','Color',colors(ii,:));
end

%% linearity of surround...
figure(2);
subplot(131); hold on;
limUp = max([tempResponse.LinearCenterPlusNonlinearSurround,tempResponse.LNCenterSurround]);
plot([0 limUp],[0 limUp],'k--')
plot(tempResponse.LinearCenterPlusNonlinearSurround,tempResponse.LNCenterSurround,...
    'Marker','o','LineStyle','none','Color',colors(ii,:));

subplot(132); hold on;
limUp = max([tempResponse.SubunitCenterPlusSubunitSurround,tempResponse.SubunitCenterPlusLinearSurround]);
plot([0 limUp],[0 limUp],'k--')
plot(tempResponse.SubunitCenterPlusSubunitSurround,tempResponse.SubunitCenterPlusLinearSurround,...
    'Marker','o','LineStyle','none','Color',colors(ii,:));


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
for ss = 1:length(spotSize) %get responses to each spot
    currentRadius = spotSize(ss)/2;
    spotBinary = double(sqrt((rr-(patchSize/2)).^2+(cc-(patchSize/2)).^2)<=currentRadius);
    responseStructure = RFmodel.getResponse(spotBinary);
    
    response.LNCenterSurround(ss) = responseStructure.CenterSurround.LN;
    response.SubunitCenterPlusLinearSurround(ss) = responseStructure.CenterSurround.SharedNonlinearity;
    response.SubunitCenterPlusSubunitSurround(ss) = responseStructure.CenterSurround.LNLNSamePolarity;
    response.LinearCenterPlusNonlinearSurround(ss) = responseStructure.CenterSurround.LinearCenterPlusNonlinearSurround;
end
spotSize = spotSize .* MicronsPerPixel;
figure(2); clf; hold on
plot(spotSize,response.LNCenterSurround ./ max(response.LNCenterSurround),'b-')
plot(spotSize,response.SubunitCenterPlusLinearSurround ./ max(response.SubunitCenterPlusLinearSurround),'r-o')
plot(spotSize,response.SubunitCenterPlusSubunitSurround ./ max(response.SubunitCenterPlusSubunitSurround),'k-')
plot(spotSize,response.LinearCenterPlusNonlinearSurround ./ max(response.LinearCenterPlusNonlinearSurround),'g-')

%%
% gratings...
barSizes = [5 10 20 30 40 50 60 70 80 90 100,...
    120 140 160 180 200 220 240 260 280 300 320];
response.LNCenterSurround = [];
response.SubunitCenterPlusLinearSurround = [];
response.SubunitCenterPlusSubunitSurround = [];
response.LinearCenterPlusNonlinearSurround = [];
for bb = 1:length(barSizes);
    barSize = barSizes(bb); %microns
    sf = 1 / (2*barSize/MicronsPerPixel);
    tempXX = sin(2*pi.*sf*(-patchSize/2 + 1: patchSize/2));
    newStim = repmat(tempXX,patchSize,1);
    responseStructure = RFmodel.getResponse(newStim);
    
    response.LNCenterSurround(bb) = responseStructure.CenterSurround.LN;
    response.SubunitCenterPlusLinearSurround(bb) = responseStructure.CenterSurround.SharedNonlinearity;
    response.SubunitCenterPlusSubunitSurround(bb) = responseStructure.CenterSurround.LNLNSamePolarity;
    response.LinearCenterPlusNonlinearSurround(bb) = responseStructure.CenterSurround.LinearCenterPlusNonlinearSurround;
end

figure(2); clf; hold on
plot(barSizes,response.LNCenterSurround ./ max(response.LNCenterSurround),'b-')
plot(barSizes,response.SubunitCenterPlusLinearSurround ./ max(response.SubunitCenterPlusLinearSurround),'r-o')
plot(barSizes,response.SubunitCenterPlusSubunitSurround ./ max(response.SubunitCenterPlusSubunitSurround),'k-')
plot(barSizes,response.LinearCenterPlusNonlinearSurround ./ max(response.LinearCenterPlusNonlinearSurround),'g-')
