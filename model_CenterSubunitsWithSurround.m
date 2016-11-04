scaleFactor = 3.3; %um per arcmin

barWidths = round((10:10:300) ./ scaleFactor);
gratingContrast = 0.9;


% RF properties
FilterSize = round(600 ./ scaleFactor); % microns -> arcmin

%stdevs... 
CenterSigma = round(40 ./ scaleFactor);
SubunitCenterSigma = round((CenterSigma/6));
SubunitSurroundSigma = round(CenterSigma);
SubunitCenterWeight = 0.6;

  
% subunit locations - square grid
TempFilter = zeros(FilterSize, FilterSize);
SubunitLocations = find(rem([1:FilterSize], 2*SubunitCenterSigma) == 0);
for x = 1:length(SubunitLocations)
    TempFilter(SubunitLocations(x), SubunitLocations) = 1;
end
SubunitIndices = find(TempFilter > 0);

% center & subunit filters
for x = 1:FilterSize
    for y = 1:FilterSize
        SubunitCenter(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (SubunitCenterSigma^2)));
        SubunitSurround(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (SubunitSurroundSigma^2)));
        RFCenter(x,y) = exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) / (2 * (CenterSigma^2)));
    end
end
subunitWeightings = RFCenter(SubunitIndices);


% normalize each component
subunitWeightings = subunitWeightings / sum(subunitWeightings);
SubunitCenter = SubunitCenter / sum(SubunitCenter(:));
SubunitSurround = SubunitSurround / sum(SubunitSurround(:));

LinearSubunitFilter = SubunitCenterWeight * SubunitCenter - (1 - SubunitCenterWeight) * SubunitSurround;

%get RF activations for each grating
RFoutput_ln = [];
RFoutput_subunit = [];

for bb = 1:length(barWidths)
    tempStim_right = sin(2*(1/(2*barWidths(bb)))*pi.*(0:FilterSize/2-1));
    tempStim_left = sin(2*(1/(2*barWidths(bb)))*pi.*(1:FilterSize/2));
    
    tempStim = [-fliplr(tempStim_left) tempStim_right];
    
    tempStim(tempStim>0) = gratingContrast; tempStim(tempStim<=0) = -gratingContrast; %contrast
    stimulusImage = repmat(tempStim,FilterSize,1);

    % convolve patch with subunit filter
    ImagePatch = conv2(stimulusImage, LinearSubunitFilter, 'same');

    % activation of each subunit
    subunitActivations = ImagePatch(SubunitIndices);

    %nonlinear subunits - rect linear as subunit nlinearity
    subunitOutputs = subunitActivations;
    subunitOutputs(subunitOutputs<0) = 0;

    % Linear center
    RFoutput_ln(bb) = sum(subunitActivations .* subunitWeightings);
    %subunit...
    RFoutput_subunit(bb) = sum(subunitOutputs.* subunitWeightings);

end
%rectify LN output...
LnModelResponse = RFoutput_ln;
LnModelResponse(LnModelResponse < 0) = 0;

SubunitModelResponse = RFoutput_subunit;

figure(1); clf;
plot(barWidths .* scaleFactor, LnModelResponse ,'k')
hold on;
plot(barWidths .* scaleFactor, SubunitModelResponse ,'r')



