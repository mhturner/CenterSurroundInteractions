function movieMatrix = getDovesMovie(stimulusIndex,frameSize)


resourcesDir = '~/Documents/MATLAB/turner-package/resources/';
imageSet = 'NaturalImages';
stimSet = 'dovesFEMstims_20160422.mat';
load([resourcesDir, stimSet]);

imageName = FEMdata(stimulusIndex).ImageName;
fileId=fopen([resourcesDir, imageSet, '/', imageName],'rb','ieee-be');
img = fread(fileId, [1536,1024], 'uint16');
img = double(img);
img = (img./max(img(:))); %rescale s.t. brightest point is maximum monitor level

noFrames = length(FEMdata(stimulusIndex).eyeX);
movieMatrix = zeros(frameSize,frameSize,noFrames);
for ff = 1:noFrames
    centerX = FEMdata(stimulusIndex).eyeX(ff);
    centerY = FEMdata(stimulusIndex).eyeY(ff);
    movieMatrix(:,:,ff) = img(round(centerX - frameSize/2 + 1) : round(centerX + frameSize/2),...
        round(centerY - frameSize/2 + 1) : round(centerY + frameSize/2));
end

end