classdef SpatioTemporalRFModel < handle
    
    properties
        MicronsPerPixel = 6.6; %microns
        
        %spatial RF components
        subunitSurroundWeight = 0.72; %As a fraction of subunit center weight
        filterSize_um = 700;
        subunitSigma_um = 10;
        subunitSurroundSigma_um = 150;
        centerSigma_um = 40;
        
        %temporal RF components
        surroundFilterTimeShift = 0; %msec
    end
    
    properties
        SubunitFilter
        SubunitSurroundFilter
        SubunitIndices
        SubunitWeightings
        
        CenterTemporalFilter
        SurroundTemporalFilter
    end

    methods
        
        function response = getResponse(obj,stimulusMatrix,stimulusWave)
            if isempty(obj.SubunitFilter)
               error('Initialize model with makeRfComponents')
            end
            if ~isequal(size(stimulusMatrix),size(obj.SubunitFilter))
                error(['stimulusMatrix must be size: ', ...
                    num2str(size(obj.SubunitFilter,1)), ' by ',num2str(size(obj.SubunitFilter,2))])
            end
            
            %activation in space:
            convolved_SubunitCenter = conv2(stimulusMatrix, obj.SubunitFilter, 'same');
            subunitCenterActivations = convolved_SubunitCenter(obj.SubunitIndices);

            convolved_SubunitSurround = conv2(stimulusMatrix, obj.SubunitSurroundFilter, 'same');
            subunitSuroundActivations = convolved_SubunitSurround(obj.SubunitIndices);

            %temporal activation of each subunit center & surround
            %rows = subunit, cols = time
            centerFilteredStim = conv(stimulusWave,obj.CenterTemporalFilter);
            centerFilteredStim = centerFilteredStim(1:length(stimulusWave));
            
            surroundFilteredStim = conv(stimulusWave,obj.SurroundTemporalFilter);
            surroundFilteredStim = surroundFilteredStim(1:length(stimulusWave));
            
            subunitVoltage = subunitCenterActivations * centerFilteredStim -...
                (obj.subunitSurroundWeight .* subunitSuroundActivations) * surroundFilteredStim;

            %hit each subunit with rectifying synapse
            subunitVoltage(subunitVoltage < 0) = 0;
            subunitOutputs = subunitVoltage;
            response = sum(subunitOutputs .* obj.SubunitWeightings); %sum over subunits
        end
        
        function makeRfComponents(obj)
            resources_dir = '~/Dropbox/RiekeLab/Analysis/MATLAB/CenterSurroundInteractions/resources/';
            % % % % % % % % MAKE SPATIAL RF COMPONENTS % % % % % % % % % % % % % % % % 
            %convert to pixels:
            filterSize = obj.Micron2Pixel(obj.filterSize_um);
            subunitSigma = obj.Micron2Pixel(obj.subunitSigma_um);
            subunitSurroundSigma = obj.Micron2Pixel(obj.subunitSurroundSigma_um);
            centerSigma = obj.Micron2Pixel(obj.centerSigma_um);

            TempFilter = zeros(filterSize, filterSize);
            SubunitLocations = find(rem((1:filterSize), 2*subunitSigma) == 0);
            for x = 1:length(SubunitLocations)
                TempFilter(SubunitLocations(x), SubunitLocations) = 1;
            end
            obj.SubunitIndices = find(TempFilter > 0);

            % Filters for subunit with surround, center, surround
            for x = 1:filterSize
                for y = 1:filterSize
                    obj.SubunitFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (subunitSigma^2)));
                    obj.SubunitSurroundFilter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (subunitSurroundSigma^2)));
                    RFCenter(x,y) = exp(-((x - filterSize/2).^2 + (y - filterSize/2).^2) / (2 * (centerSigma^2)));
                end
            end
            % normalize components
            obj.SubunitFilter = obj.SubunitFilter / sum(obj.SubunitFilter(:));
            obj.SubunitSurroundFilter = obj.SubunitSurroundFilter / sum(obj.SubunitSurroundFilter(:));
            RFCenter = RFCenter / sum(RFCenter(:));
            %get weighting of each subunit output by center
            obj.SubunitWeightings = RFCenter(obj.SubunitIndices);
            obj.SubunitWeightings = obj.SubunitWeightings ./ sum(obj.SubunitWeightings);
            
            % % % % % % % % MAKE TEMPORAL RF COMPONENTS % % % % % % % % % % % % % % % 
            load([resources_dir, 'OffParasolExcitatoryFilters.mat'])
            obj.CenterTemporalFilter = mean(filters.center);
            obj.CenterTemporalFilter = decimate(obj.CenterTemporalFilter,10); %resample to msec from 10 Khz
            obj.CenterTemporalFilter = -obj.CenterTemporalFilter ./ max(abs(obj.CenterTemporalFilter));
            
            obj.SurroundTemporalFilter = mean(filters.surround);
            obj.SurroundTemporalFilter = decimate(obj.SurroundTemporalFilter,10); %resample to msec from 10 Khz
            obj.SurroundTemporalFilter = obj.SurroundTemporalFilter ./ max(abs(obj.SurroundTemporalFilter));
            
            if obj.surroundFilterTimeShift > 0 %shift surround later in time
                obj.SurroundTemporalFilter = cat(2,zeros(1,obj.surroundFilterTimeShift),obj.SurroundTemporalFilter);
                obj.SurroundTemporalFilter((length(obj.CenterTemporalFilter) + 1) : end) = [];
            elseif obj.surroundFilterTimeShift < 0 %shift surround earlier in time
                obj.SurroundTemporalFilter(1:abs(obj.surroundFilterTimeShift)) = [];
                obj.SurroundTemporalFilter(end:length(obj.CenterTemporalFilter)) = 0;
            else
                %no shift
            end
        end

    end
    methods (Access = private)
        function res = Micron2Pixel(obj,microns)
            res = round(microns / obj.MicronsPerPixel);
        end
    end
    
end