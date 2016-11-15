clear all; close all; clc;
load('NaturalImageFlashLibrary_101716.mat')
imageNames = fieldnames(imageData);
resourcesDir = '~/Documents/MATLAB/turner-package/resources';
patchSize = 100; %pixels
MicronsPerPixel = 6.6;
patchesPerImage = 5e2;
noBins = 8;
noImages = 5; %up to 20

RFmodel = ThreeLayerReceptiveFieldModel;
RFmodel.MicronsPerPixel = MicronsPerPixel;
RFmodel.SubunitSurroundWeight = 0.9;
RFmodel.SurroundSubunitSigma = 12;
RFmodel.SubunitSurroundSamplingSigma = 150;
RFmodel.makeRfComponents(patchSize);


ImageX = 1536; ImageY = 1024;
distanceBins = linspace(0,500,noBins+1);
% distanceBins = linspace(0,sqrt((ImageY-patchSize)^2 + (ImageX-patchSize)^2),noBins+1);
binCenters = distanceBins(1:end-1) + diff(distanceBins)./2;
binCenters = binCenters .* MicronsPerPixel; %um

colors = pmkmp(20);
corrFieldNames = {'imageContrast','LNcenter','subunitCenter','LNCenterSurround',...
    'SubunitCenterPlusLinearSurround','SubunitCenterPlusSubunitSurround'};
correlationFunction = struct;
for cc = 1:length(corrFieldNames)
    correlationFunction.(corrFieldNames{cc}) = zeros(noImages, noBins);
end

for ii = 1:noImages
    disp(num2str(ii))
    fileId=fopen([resourcesDir, '/VHsubsample_20160105', '/', imageNames{ii},'.iml'],'rb','ieee-be');
    img = fread(fileId, [1536,1024], 'uint16');
    contrastImg = img - mean(img(:));
    contrastImg = contrastImg ./ max(contrastImg(:));
    
    location = zeros(patchesPerImage,2);
    tempResponse = struct;
    for cc = 1:length(corrFieldNames)
        tempResponse.(corrFieldNames{cc}) = zeros(patchesPerImage,1);
    end
    for rr = 1:patchesPerImage
         % choose location
        x = round(patchSize/2 + (ImageX - patchSize)*rand);
        y = round(patchSize/2 + (ImageY - patchSize)*rand);
        location(rr,:) = [x, y];

        % get patch
        ContrastPatch = contrastImg(x-patchSize/2+1:x+patchSize/2,y-patchSize/2+1:y+patchSize/2);
        
        % measure spatial contrast
        tempResponse.imageContrast(rr) = var(ContrastPatch(:));
        % get RF model responses
        responseStructure = RFmodel.getResponse(ContrastPatch);
        tempResponse.LNcenter(rr) = responseStructure.CenterOnly.LN;
        tempResponse.subunitCenter(rr) = responseStructure.CenterOnly.NonlinearSubunits;
        
        tempResponse.LNCenterSurround(rr) = responseStructure.CenterSurround.LN;
        tempResponse.SubunitCenterPlusLinearSurround(rr) = responseStructure.CenterSurround.SharedNonlinearity;
        tempResponse.SubunitCenterPlusSubunitSurround(rr) = responseStructure.CenterSurround.LNLNSamePolarity;
    end
    distanceMatrix = squareform(pdist(location));
    
    for cc = 1:length(corrFieldNames)
        tempCorrs = zeros(1,noBins);
        for dd = 1:noBins
            reference = []; afar = [];
            for pp = 1:patchesPerImage
                distances = distanceMatrix(pp,:);
                inds = find(distances > distanceBins(dd) & distances < distanceBins(dd+1));
                afar = cat(1,afar,tempResponse.(corrFieldNames{cc})(inds));
                reference = cat(1,reference,tempResponse.(corrFieldNames{cc})(pp).*ones(length(inds),1));
            end
            tempCorrs(dd) = corr(reference,afar);
        end
        correlationFunction.(corrFieldNames{cc})(ii,:) = tempCorrs;
    end
    
    
end
%%

corrFieldNames = {'imageContrast','LNcenter','subunitCenter','LNCenterSurround',...
    'SubunitCenterPlusLinearSurround','SubunitCenterPlusSubunitSurround'};

figure(3); clf; hold on;
xlabel('Distance (um)'); ylabel('Correlation')

%image contrast
meanCorr = mean(correlationFunction.imageContrast,1);
h1 = plot(binCenters,meanCorr,'k:');
plot([0 binCenters(end)],[0 0],'k--')
% for bb = 1:length(binCenters)
%     err = std(correlationFunction.imageContrast(:,bb));
%     plot([binCenters(bb) binCenters(bb)],...
%         meanCorr(bb) + [err -err],'k-') 
% end

%LN center
meanCorr = mean(correlationFunction.LNcenter,1);
h2 = plot(binCenters,meanCorr,'b:o');
% for bb = 1:length(binCenters)
%     err = std(correlationFunction.LNcenter(:,bb));
%     plot([binCenters(bb) binCenters(bb)],...
%         meanCorr(bb) + [err -err],'b-') 
% end

%Sub center
meanCorr = mean(correlationFunction.subunitCenter,1);
h3 = plot(binCenters,meanCorr,'r:o');
% for bb = 1:length(binCenters)
%     err = std(correlationFunction.subunitCenter(:,bb));
%     plot([binCenters(bb) binCenters(bb)],...
%         meanCorr(bb) + [err -err],'r-') 
% end

%LN CenterSurround
meanCorr = mean(correlationFunction.LNCenterSurround,1);
h4 = plot(binCenters,meanCorr,'b-o');
% for bb = 1:length(binCenters)
%     err = std(correlationFunction.LNCenterSurround(:,bb));
%     plot([binCenters(bb) binCenters(bb)],...
%         meanCorr(bb) + [err -err],'b--') 
% end

%SubunitCenterPlusLinearSurround
meanCorr = mean(correlationFunction.SubunitCenterPlusLinearSurround,1);
h5 = plot(binCenters,meanCorr,'r-o');
% for bb = 1:length(binCenters)
%     err = std(correlationFunction.SubunitCenterPlusLinearSurround(:,bb));
%     plot([binCenters(bb) binCenters(bb)],...
%         meanCorr(bb) + [err -err],'r--') 
% end

%SubunitCenterPlusSubunitSurround
meanCorr = mean(correlationFunction.SubunitCenterPlusSubunitSurround,1);
h6 = plot(binCenters,meanCorr,'g-o');
% for bb = 1:length(binCenters)
%     err = std(correlationFunction.SubunitCenterPlusSubunitSurround(:,bb));
%     plot([binCenters(bb) binCenters(bb)],...
%         meanCorr(bb) + [err -err],'g--') 
% end

legend([h1 h2 h3 h4 h5 h6],'Image: spatial contrast','LN center','Subunit center',...
    'LN center-surround','Subunit center + linear surround','Subunit center + subunit surround')
    




%%
figure(1); clf;



meanPhi = [];
i0 = [0:10:1000];
for ii = 1:length(i0)
    phi = log(img(:)./i0(ii));
    meanPhi(ii) = mean(phi);
end
    
plot(i0,meanPhi,'k')
hold on;
plot([i0(1) i0(end)],[0 0],'k--')

%%


[N,edges] = histcounts(phi);
p = N ./ sum(N);
binCenters = edges(1:end-1) + diff(edges);

figure(2); clf;
semilogy(binCenters, p,'k')

%%

figure(3); clf;
subplot(221)
imagesc(img'); colormap(gray); axis image;
[Gmag,~] = imgradient(img);
subplot(222)
imagesc(Gmag'); colormap(gray); axis image;

phi = log(img./827);

subplot(223)
imagesc(phi'); colormap(gray); axis image;
[Gmag,~] = imgradient(phi);
subplot(224)
imagesc(Gmag'); colormap(gray); axis image;


[N,edges] = histcounts(Gmag(:));
p = N ./ sum(N);
binCenters = edges(1:end-1) + diff(edges);
figure(4); clf; 
semilogy(binCenters, p,'k')

vnImg = [];
%variance-normalized image:
N = 4;
for yy = (N/2):size(img,1)-(N/2)
    for xx = (N/2):size(img,2)-(N/2)
        patch = img((yy-N/2+1):(yy+N/2),(xx-N/2+1):(xx+N/2));
        vnImg(xx,yy) = (img(yy,xx) - mean(patch(:))) / std(patch(:));
    end
end

phi = log(vnImg./827);
[Gmag,~] = imgradient(phi);
[N,edges] = histcounts(Gmag(:));
p = N ./ sum(N);
binCenters = edges(1:end-1) + diff(edges);

hold on;
semilogy(binCenters, p,'b')