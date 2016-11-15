classdef ThreeLayerReceptiveFieldModel < handle
    
    properties
        SurroundSubunitSigma = 24; % (microns) gaussian subunits, samples pixels (~horizontal cell subunits)
        
        CenterSubunitSigma = 12; %gaussian subunits, samples pixels (~bipolar cell center)
        SubunitSurroundSamplingSigma = 120; %sigma, samples surround units (~bipolar cell surround)
        SubunitSurroundWeight = 0.8; %As a fraction of subunit center weight

        CenterSamplingSigma = 40; %samples layer 2 subunits
        
        MicronsPerPixel = 3.3; %microns

    end
    
    properties
        convolutionFilter
        samplingFilter
        sampleWeights
        
        
    end
    
    properties (Hidden)
        SubunitActivation
        SubunitIndices
    end
    
    methods
        
        function responseStructure = getResponse(obj,stimulusMatrix)
            if isempty(obj.convolutionFilter)
               error('Initialize model with makeRfComponents')
            end

            obj.getSubunitActivations(stimulusMatrix);

            %Linear component responses:
            LinearCenterActivation = sum(obj.SubunitActivation.CenterSubunits .* obj.sampleWeights.Center);
            LinearSurroundActivation = sum(obj.SubunitActivation.SurroundSubunits .* obj.sampleWeights.SubunitSurround);
            
            responseStructure.component.LinearCenterActivation = LinearCenterActivation;
            responseStructure.component.LinearSurroundActivation = LinearSurroundActivation;
            
            responseStructure.component.NonlinearCenterActivation = max(0,LinearCenterActivation);
            responseStructure.component.NonlinearSurroundActivation_samePolarity = max(0,LinearSurroundActivation);

            % % % % % % % RF CENTER ONLY % % % % % % % % 
            % 1) LN Center:
            responseStructure.CenterOnly.LN = max(0,LinearCenterActivation);
            % 2) Nonlinear subunits in the center:
            subunitOutputs = obj.SubunitActivation.CenterSubunits;
            subunitOutputs(subunitOutputs < 0) = 0; %rectify subunit outputs
            responseStructure.CenterOnly.NonlinearSubunits =...
                sum(subunitOutputs.* obj.sampleWeights.Center);
            
            
            % % % % % % % ADD RF SURROUND % % % % % % % % 
            % 3) R = N(Linear surround + linear center)
            responseStructure.CenterSurround.LN = ...
                max(0, LinearCenterActivation - LinearSurroundActivation);
            % 4) R = N(Center) + N(Surround)
            subunitOutputs = obj.SubunitActivation.SurroundSubunits;
            subunitOutputs(subunitOutputs < 0) = 0; %rectify surround subunit outputs
            responseStructure.CenterSurround.NonlinearCenterPlusSurround = ...
                responseStructure.CenterOnly.NonlinearSubunits - ...
                sum(subunitOutputs.* obj.sampleWeights.SubunitSurround);
            
            % 5) R = N(Center + surround) - i.e. subunitsHaveSurround
            SurroundInputToEachSubunit = ...
                obj.sampleWeights.SubunitWeightMatrix * obj.SubunitActivation.SurroundSubunits;
            subunitActivation = obj.SubunitActivation.CenterSubunits - SurroundInputToEachSubunit;
            subunitOutputs = subunitActivation;
            subunitOutputs(subunitOutputs < 0) = 0; %rectify subunit outputs
            responseStructure.CenterSurround.SharedNonlinearity =...
                sum(subunitOutputs.* obj.sampleWeights.Center);
            
            % 6) Surround subunit output is rectified, then combined with
            % center subunits before passing thru a shared nlinearity
            % (LNLN)
            SurroundInputToEachSubunit = ...
                max(0,obj.sampleWeights.SubunitWeightMatrix * obj.SubunitActivation.SurroundSubunits);
            subunitActivation = obj.SubunitActivation.CenterSubunits - SurroundInputToEachSubunit;
            subunitOutputs = subunitActivation;
            subunitOutputs(subunitOutputs < 0) = 0; %rectify subunit outputs
            responseStructure.CenterSurround.LNLNSamePolarity =...
                sum(subunitOutputs.* obj.sampleWeights.Center);
            
            % 7) R = N(linear center - nonlinear surround)
            SurroundInputToEachSubunit = ...
                max(0,obj.sampleWeights.SubunitWeightMatrix * obj.SubunitActivation.SurroundSubunits);
            TotalSurroundInput = sum(SurroundInputToEachSubunit .* obj.sampleWeights.Center);
            responseStructure.CenterSurround.LinearCenterPlusNonlinearSurround = ...
                max(0, LinearCenterActivation - TotalSurroundInput);
            
            % DIVISIVE models:
            additiveConstant = 0.01;
            % linCenter & surround
            linCenter = max(0, LinearCenterActivation);
            linSurround = max(0,LinearSurroundActivation);
            % nlCenter
            nlCenter = responseStructure.CenterOnly.NonlinearSubunits;
            % nlSurround
            SurroundInputToEachSubunit = ...
                max(0,obj.sampleWeights.SubunitWeightMatrix * obj.SubunitActivation.SurroundSubunits);
            nlSurround = sum(SurroundInputToEachSubunit .* obj.sampleWeights.Center);
            
            responseStructure.DivisiveSurround.linC_linS = linCenter / (linSurround + additiveConstant);
            responseStructure.DivisiveSurround.linC_nlS = linCenter / (nlSurround + additiveConstant);
            responseStructure.DivisiveSurround.nlC_linS = nlCenter / (linSurround + additiveConstant);
            responseStructure.DivisiveSurround.nlC_nlS = nlCenter / (nlSurround + additiveConstant);
            
        end
        
        function makeRfComponents(obj,FilterSize)
            % Filter size in pixels
            % make subunit location grids and subunit indices to sample
            % pixels
            %    1) surround (horizontal) subunits:
            SurroundSubunitLocations = find(rem(1:FilterSize,...
               2 * obj.Micron2Pixel(obj.SurroundSubunitSigma)) == 0);
            TempFilter = zeros(FilterSize, FilterSize);
            for x = 1:length(SurroundSubunitLocations)
                TempFilter(SurroundSubunitLocations(x), SurroundSubunitLocations) = 1;
            end
            obj.SubunitIndices.Surround = find(TempFilter > 0);
            %    2) center (bipolar) subunits:
            CenterSubunitLocations = find(rem(1:FilterSize,...
               2 * obj.Micron2Pixel(obj.CenterSubunitSigma)) == 0);
            TempFilter = zeros(FilterSize, FilterSize);
            for x = 1:length(CenterSubunitLocations)
                TempFilter(CenterSubunitLocations(x), CenterSubunitLocations) = 1;
            end
            obj.SubunitIndices.Center = find(TempFilter > 0);
            
            % Make filters:
            for x = 1:FilterSize
                for y = 1:FilterSize
                    % Two convolution filters:
                    obj.convolutionFilter.SurroundSubunit(x,y) = ...
                        exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) /...
                        (2 * (obj.Micron2Pixel(obj.SurroundSubunitSigma)^2)));
                    obj.convolutionFilter.CenterSubunit(x,y) = ...
                        exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) /...
                        (2 * (obj.Micron2Pixel(obj.CenterSubunitSigma)^2)));
                    
                    % Two sampling filters:
                    obj.samplingFilter.SubunitSurround(x,y) = ...
                        exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) /...
                        (2 * (obj.Micron2Pixel(obj.SubunitSurroundSamplingSigma)^2)));
                    obj.samplingFilter.Center(x,y) = ...
                        exp(-((x - FilterSize/2).^2 + (y - FilterSize/2).^2) /...
                        (2 * (obj.Micron2Pixel(obj.CenterSamplingSigma)^2)));
                end
            end
            % get arrays of sample weights-
            obj.sampleWeights.Center = obj.samplingFilter.Center(obj.SubunitIndices.Center);
            obj.sampleWeights.SubunitSurround = obj.samplingFilter.SubunitSurround(obj.SubunitIndices.Surround);

            % Get the subunit weight matrix. This is size (number of center
            % subunits) by (number of surround subunits). Each row defines the
            % connection weight from each surround subunit to that center
            % subunit. Each row sums to obj.SubunitSurroundWeight
            obj.sampleWeights.SubunitWeightMatrix = ...
                zeros(length(obj.SubunitIndices.Center),length(obj.SubunitIndices.Surround));
            tempMat = zeros(FilterSize, FilterSize);
            for cc = 1:length(obj.SubunitIndices.Center)
                currentMat = tempMat;
                currentMat(obj.SubunitIndices.Center(cc)) = 1;
                currentMat = conv2(currentMat,obj.samplingFilter.SubunitSurround,'same');
                tempRow = currentMat(obj.SubunitIndices.Surround);
                %normalize row to obj.SubunitSurroundWeight
                tempRow = obj.SubunitSurroundWeight .* tempRow ./ sum(tempRow(:));
                obj.sampleWeights.SubunitWeightMatrix(cc,:) = tempRow;
            end

            % normalize components:
            %   1) convolution filters
            obj.convolutionFilter.SurroundSubunit = ...
                obj.convolutionFilter.SurroundSubunit / sum(obj.convolutionFilter.SurroundSubunit(:)); 
            obj.convolutionFilter.CenterSubunit = ...
                obj.convolutionFilter.CenterSubunit / sum(obj.convolutionFilter.CenterSubunit(:));
            %   2) Sample weights
            %sums to obj.SubunitSurroundWeight:
            obj.sampleWeights.SubunitSurround = ...
                obj.SubunitSurroundWeight .* obj.sampleWeights.SubunitSurround /...
                sum(obj.sampleWeights.SubunitSurround(:)); 
            %sums to one:
            obj.sampleWeights.Center = ...
                obj.sampleWeights.Center / sum(obj.sampleWeights.Center(:));
        end
        
        function getSubunitActivations(obj,stimulusMatrix)
            % Convolution steps - independent of sampling architecture:
            % 1) Surround subunit activations:
            surroundStimulus_convolved = ...
                conv2(stimulusMatrix, obj.convolutionFilter.SurroundSubunit, 'same');
            % activation of each subunit
            obj.SubunitActivation.SurroundSubunits = surroundStimulus_convolved(obj.SubunitIndices.Surround);
            
            % 1) Center subunit activations:
            centerStimulus_convolved = ...
                conv2(stimulusMatrix, obj.convolutionFilter.CenterSubunit, 'same');
            % activation of each subunit
            obj.SubunitActivation.CenterSubunits = centerStimulus_convolved(obj.SubunitIndices.Center);
        end

        


    end
    methods (Access = private)
        function res = Micron2Pixel(obj,microns)
            res = round(microns / obj.MicronsPerPixel);
        end
    end
    
end