function doCSLNAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'bins2D',6^2,@isnumeric);
    addParameter(ip,'bins1D',20,@isnumeric); 
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'convertToConductance',true,@islogical);
    addParameter(ip,'fitWithEquallyPopulatedBins',true,@islogical);
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    bins2D = ip.Results.bins2D;
    bins1D = ip.Results.bins1D;
    exportFigs = ip.Results.exportFigs;
    convertToConductance = ip.Results.convertToConductance;
    fitWithEquallyPopulatedBins = ip.Results.fitWithEquallyPopulatedBins;
    
    figColors = pmkmp(8);
    egEpochToPull = 3; %for example traces

    figure; clf; fig1=gca; %Linear filters
    set(fig1,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig1,'XLabel'),'String','Time (s)')
    set(get(fig1,'YLabel'),'String','')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig2=gca; %Nonlinearities
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Linear prediction (nS)')
    set(get(fig2,'YLabel'),'String','Measured (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig3=gca; %Independent nonlinearities
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Linear prediction (nS)')
    set(get(fig3,'YLabel'),'String','Measured (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig4=gca; %Shared nonlinearity
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','Linear prediction (nS)')
    set(get(fig4,'YLabel'),'String','Measured (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig5=gca; %population R-squared results
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','R^2 independent nonlinearities')
    set(get(fig5,'YLabel'),'String','R^2 shared nonlinearity')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig6=gca; %slice of center, modulate surround
    set(fig6,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig6,'XLabel'),'String','Surround activation')
    set(get(fig6,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig7=gca; %slice of surround, modulate center
    set(fig7,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig7,'XLabel'),'String','Center activation')
    set(get(fig7,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig8=gca; %center +/- surround overlay
    set(fig8,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig8,'XLabel'),'String','Center +/- surround')
    set(get(fig8,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig9=gca; %surround +/- center overlay
    set(fig9,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig9,'XLabel'),'String','Surround +/- center')
    set(get(fig9,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    % Individual traces for example cell:
    figure; clf; fig11=gca; %center stim
    set(fig11,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig11,'XLabel'),'String','Time (s)')
    set(get(fig11,'YLabel'),'String','Contrast')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig12=gca; %surround stim
    set(fig12,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig12,'XLabel'),'String','Time (s)')
    set(get(fig12,'YLabel'),'String','Contrast')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig13=gca; %center response
    set(fig13,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig13,'XLabel'),'String','Time (s)')
    set(get(fig13,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig14=gca; %surround response
    set(fig14,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig14,'XLabel'),'String','Time (s)')
    set(get(fig14,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig15=gca; %center-surround response
    set(fig15,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig15,'XLabel'),'String','Time (s)')
    set(get(fig15,'YLabel'),'String','Response (nS)')
    set(gcf, 'WindowStyle', 'docked')
    
    figure; clf; fig16=gca; %improvement in R2 vs surround weight
    set(fig16,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig16,'XLabel'),'String','Relative surround weight')
    set(get(fig16,'YLabel'),'String','Improvement over independent model')
    set(gcf, 'WindowStyle', 'docked')

    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(useRandomSeed)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    
    filters.center = [];
    filters.surround = [];
    
    rSquaredValues.joint = [];
    rSquaredValues.indep = [];
    rSquaredValues.shared = [];
    ONcellInds = [];
    OFFcellInds = [];
    for pp = 1:length(populationNodes)
        recNode = populationNodes{pp};
        currentNode = [];
        repeatedNode = [];
        for cc = 1:recNode.children.length
           if recNode.children(cc).splitValue == 1 %use only random seed for now
               currentNode = recNode.children(cc);
           elseif recNode.children(cc).splitValue == 0
               repeatedNode = recNode.children(cc);
           end
        end
        cellInfo = getCellInfoFromEpochList(currentNode.epochList);
        recType = getRecordingTypeFromEpochList(currentNode.epochList);
        if (convertToConductance)
            recType = [recType, ', conductance']; %#ok<AGROW>
        end
        
        if strcmp(cellInfo.cellType,'ONparasol')
            ONcellInds = cat(2,ONcellInds,pp);
        elseif strcmp(cellInfo.cellType,'OFFparasol')
            OFFcellInds = cat(2,OFFcellInds,pp);
        end
        
% % % % % % % % GET LINEAR FILTERS AND GENERATOR SIGNALS % % % % % % % % % % % % % % % % 
        % center:
        center = getLinearFilterAndPrediction(currentNode.childBySplitValue('Center').epochList,recType,...
            'seedName','centerNoiseSeed','numberOfBins',bins1D);
        % surround:
        surround = getLinearFilterAndPrediction(currentNode.childBySplitValue('Surround').epochList,recType,...
            'seedName','surroundNoiseSeed','numberOfBins',bins1D);
        % center + surround: just to test against measured response
        centerSurround = getLinearFilterAndPrediction(currentNode.childBySplitValue('Center-Surround').epochList,recType,...
            'seedName','centerNoiseSeed','numberOfBins',bins1D);
        
        filters.center(pp,:) = center.LinearFilter;
        filters.surround(pp,:) = surround.LinearFilter;

        if currentNode.custom.get('isExample')
                addLineToAxis([center.filterTimeVector],[center.LinearFilter],...
                    ['center',num2str(pp)],fig1,figColors(1,:),'-','none')
                addLineToAxis([surround.filterTimeVector],[surround.LinearFilter],...
                    ['surround',num2str(pp)],fig1,figColors(4,:),'-','none')
                
                addLineToAxis(center.nonlinearity.binMean,center.nonlinearity.respMean,...
                    ['centerDataMean',num2str(pp)],fig2,figColors(1,:),'none','o')
                addLineToAxis(center.nonlinearity.binMean,...
                    center.nonlinearity.respMean + center.nonlinearity.respErr,...
                    ['centerDataErrUp',num2str(pp)],fig2,figColors(1,:),'--','none')
                addLineToAxis(center.nonlinearity.binMean,...
                    center.nonlinearity.respMean - center.nonlinearity.respErr,...
                    ['centerDataErrDown',num2str(pp)],fig2,figColors(1,:),'--','none')
                addLineToAxis(center.nonlinearity.fitXX,center.nonlinearity.fitYY,...
                    ['centerFit',num2str(pp)],fig2,figColors(1,:),'-','none')
                
                addLineToAxis(surround.nonlinearity.binMean,surround.nonlinearity.respMean,...
                    ['surroundDataMean',num2str(pp)],fig2,figColors(4,:),'none','o')
                addLineToAxis(surround.nonlinearity.binMean,...
                    surround.nonlinearity.respMean + surround.nonlinearity.respErr,...
                    ['surroundDataErrUp',num2str(pp)],fig2,figColors(4,:),'--','none')
                addLineToAxis(surround.nonlinearity.binMean,...
                    surround.nonlinearity.respMean - surround.nonlinearity.respErr,...
                    ['surroundDataErrDown',num2str(pp)],fig2,figColors(4,:),'--','none')
                addLineToAxis(surround.nonlinearity.fitXX,surround.nonlinearity.fitYY,...
                    ['surroundFit',num2str(pp)],fig2,figColors(4,:),'-','none')
                
                
                cEL = sortEpochList_time(currentNode.childBySplitValue('Center').epochList);
                centerExampleEpoch = cEL.elements(egEpochToPull);
                cSeed = centerExampleEpoch.protocolSettings('centerNoiseSeed');
                centerEpochRes = getNoiseStimulusAndResponse(centerExampleEpoch,recType,'seedName','centerNoiseSeed');
                addLineToAxis(centerEpochRes.wholeTrace.timeVector(1:length(centerEpochRes.wholeTrace.stimulus)),...
                    centerEpochRes.wholeTrace.stimulus,...
                    ['stim',num2str(pp)],fig11,figColors(1,:),'-','none')
                addLineToAxis(centerEpochRes.wholeTrace.timeVector(1:length(centerEpochRes.wholeTrace.stimulus)),...
                    centerEpochRes.wholeTrace.response(1:length(centerEpochRes.wholeTrace.stimulus)),...
                    ['resp',num2str(pp)],fig13,figColors(1,:),'-','none')
                
                sEL = sortEpochList_time(currentNode.childBySplitValue('Surround').epochList);
                surroundExampleEpoch = sEL.elements(egEpochToPull);
                sSeed = surroundExampleEpoch.protocolSettings('surroundNoiseSeed');
                surroundEpochRes = getNoiseStimulusAndResponse(surroundExampleEpoch,recType,'seedName','surroundNoiseSeed');
                addLineToAxis(surroundEpochRes.wholeTrace.timeVector(1:length(surroundEpochRes.wholeTrace.stimulus)),...
                    surroundEpochRes.wholeTrace.stimulus,...
                    ['stim',num2str(pp)],fig12,figColors(4,:),'-','none')
                addLineToAxis(surroundEpochRes.wholeTrace.timeVector(1:length(surroundEpochRes.wholeTrace.stimulus)),...
                    surroundEpochRes.wholeTrace.response(1:length(surroundEpochRes.wholeTrace.stimulus)),...
                    ['resp',num2str(pp)],fig14,figColors(4,:),'-','none')
                
                sEL = sortEpochList_time(currentNode.childBySplitValue('Center-Surround').epochList);
                csExampleEpoch = sEL.elements(egEpochToPull);
                cscSeed = csExampleEpoch.protocolSettings('centerNoiseSeed');
                cssSeed = csExampleEpoch.protocolSettings('surroundNoiseSeed');
                csEpochRes = getNoiseStimulusAndResponse(csExampleEpoch,recType,'seedName','centerNoiseSeed');
                addLineToAxis(csEpochRes.wholeTrace.timeVector(1:length(csEpochRes.wholeTrace.stimulus)),...
                    csEpochRes.wholeTrace.response(1:length(csEpochRes.wholeTrace.stimulus)),...
                    ['resp',num2str(pp)],fig15,'k','-','none')
                checkCenter = ~(cscSeed == cSeed);
                checkSurround = ~(cssSeed == sSeed);
                if or(checkCenter,checkSurround)
                    warning('Example traces have different noise seeds')
                end          
        end

% % % % % % % % MODEL FITTING % % % % % % % % % % % % % % % % % % % % % % % % % 
        % Split up data into training and testing. Use repeated seed data
        % if available for testing
        if isempty(repeatedNode) %train on some random seed data. Hold one epoch out to test
            epochLen = length(centerSurround.measuredResponse) / centerSurround.n;
            epochToHold = 3; %epoch number to hold out
            testDataInds = ((epochToHold-1)*epochLen + 1):(epochToHold*epochLen);
            fitDataInds = setdiff(1:length(centerSurround.measuredResponse),testDataInds);

            trainingData.centerGS = center.generatorSignal(fitDataInds);
            trainingData.surroundGS = surround.generatorSignal(fitDataInds);
            trainingData.csMeasured = centerSurround.measuredResponse(fitDataInds);
            trainingData.cMeasured = center.measuredResponse(fitDataInds);
            trainingData.sMeasured = surround.measuredResponse(fitDataInds);
            
            testingData.centerGS = center.generatorSignal(testDataInds);
            testingData.surroundGS = surround.generatorSignal(testDataInds);
            testingData.csMeasured = centerSurround.measuredResponse(testDataInds);
        else % train on all random seed data. Test on mean repeated data.
            trainingData.centerGS = center.generatorSignal;
            trainingData.surroundGS = surround.generatorSignal;
            trainingData.csMeasured = centerSurround.measuredResponse;
            trainingData.cMeasured = center.measuredResponse;
            trainingData.sMeasured = surround.measuredResponse;
            
            % get test data...
            % center GS:
            tempCenterStruct = getLinearFilterAndPrediction(repeatedNode.childBySplitValue('Center').epochList,recType,...
                'seedName','centerNoiseSeed','numberOfBins',bins1D);
            epochLen = length(tempCenterStruct.measuredResponse) / tempCenterStruct.n;
            tempStim = tempCenterStruct.stimulus(1:epochLen);
            linearPrediction_center = conv(tempStim,center.LinearFilter);
            testingData.centerGS = linearPrediction_center(1:length(tempStim));
            % surround GS:
            tempSurroundStruct = getLinearFilterAndPrediction(repeatedNode.childBySplitValue('Surround').epochList,recType,...
                'seedName','surroundNoiseSeed','numberOfBins',bins1D);
            epochLen = length(tempSurroundStruct.measuredResponse) / tempSurroundStruct.n;
            tempStim = tempSurroundStruct.stimulus(1:epochLen);
            linearPrediction_surround = conv(tempStim,surround.LinearFilter);   
            testingData.surroundGS = linearPrediction_surround(1:length(tempStim));
            % measured c+s:
            tempCSStruct = getLinearFilterAndPrediction(repeatedNode.childBySplitValue('Center-Surround').epochList,recType,...
                'seedName','centerNoiseSeed','numberOfBins',bins1D);
            epochLen = length(tempCSStruct.measuredResponse) / tempCSStruct.n;
            responseMatrix = reshape(tempCSStruct.measuredResponse,epochLen,tempCSStruct.n)';
            testingData.csMeasured = mean(responseMatrix);
            testingData.csMeasuredVariance = var(responseMatrix);
        end
        
        %bin up and shape training data:
        if (fitWithEquallyPopulatedBins)
            %equally populated bins...
            [~,centerGS,~,centerBinID] = ...
                histcounts_equallyPopulatedBins(trainingData.centerGS,sqrt(bins2D));
            [~,surroundGS,~,surroundBinID] = ...
                histcounts_equallyPopulatedBins(trainingData.surroundGS,sqrt(bins2D));
        else
            %evenly-spaced bins...
            [~,edges,centerBinID] = histcounts(trainingData.centerGS,sqrt(bins2D));
            centerGS = edges(1:end-1) + diff(edges);
            [~,edges,surroundBinID] = histcounts(trainingData.surroundGS,sqrt(bins2D));
            surroundGS = edges(1:end-1) + diff(edges);
        end
        
        centerMean = zeros(sqrt(bins2D));
        surroundMean = zeros(sqrt(bins2D));
        responseMean = zeros(sqrt(bins2D));
        responseErr = zeros(sqrt(bins2D));
        for xx = 1:sqrt(bins2D)
            for yy = 1:sqrt(bins2D)
                jointInds = intersect(find(centerBinID == xx),find(surroundBinID == yy));
                centerMean(yy,xx) = mean(trainingData.cMeasured(jointInds));
                surroundMean(yy,xx) = mean(trainingData.sMeasured(jointInds));
                responseMean(yy,xx) = mean(trainingData.csMeasured(jointInds));
                responseErr(yy,xx) = std(trainingData.csMeasured(jointInds));
            end
        end
        %1) FIT JOINT NONLINEARITY MODEL. 7 FREE PARAMETERS
        % params = [alpha mu1 mu2 std1 std2 corr12 epsilon]
        params0 = [3*max(responseMean(:)), 400, 0, 50, 50, 0, 0];
        fitRes_joint = fitNLinearity_2D(centerGS,surroundGS,responseMean,params0);
            %2D fit surface:
            cc = linspace(min(centerGS),max(centerGS),20);
            ss = linspace(min(surroundGS),max(surroundGS),20);
            [CC,SS] = meshgrid(cc',ss');
            fitSurface = JointNLin_mvcn(CC(:)',SS(:)',fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);
            fitSurface = reshape(fitSurface,length(ss),length(cc));
            
            figure(20); clf; set(gcf, 'WindowStyle', 'docked')
            subplot(221); hold on;
            stem3(centerGS,surroundGS,responseMean)
            surf(cc,ss,fitSurface)
            xlabel('Center'); ylabel('Surroud')
            title(['Joint, ',num2str(fitRes_joint.rSquared)])
        
        % 2) FIT INDEPENDENT NONLINEARITY MODEL. 7 FREE PARAMETERS
        %params is [alphaC, betaC, gammaC,...
        %           alphaS, betaS, gammaS, epsilon]
        upGuess = 3*max(responseMean(:));
        betaGuess = center.nonlinearity.fitParams.beta;
        if pp == 7
            betaGuess = 0.1;
            upGuess = 2*max(responseMean(:));
        end
        params0 = [upGuess, betaGuess, center.nonlinearity.fitParams.gamma,...
            3*max(responseMean(:)), surround.nonlinearity.fitParams.beta, surround.nonlinearity.fitParams.gamma,...
            min([surround.nonlinearity.fitParams.epsilon, center.nonlinearity.fitParams.epsilon])];
        fitRes_indep = fitCSModel_IndependentNL(centerGS,surroundGS,responseMean,params0);
            %2D fit surface:
            fitSurface = CSModel_IndependentNL(CC(:)',SS(:)',...
                fitRes_indep.alphaC,fitRes_indep.betaC,...
                fitRes_indep.gammaC,fitRes_indep.alphaS,...
                fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);

            fitSurface = reshape(fitSurface,length(ss),length(cc));
            
            figure(20);
            subplot(222); hold off;
            stem3(centerGS,surroundGS,responseMean); hold on;
            surf(cc,ss,fitSurface)
            xlabel('Center'); ylabel('Surroud')
            title(['Indep, ',num2str(fitRes_indep.rSquared)])
            
        if currentNode.custom.get('isExample')
            xxC = min(center.generatorSignal) : max(center.generatorSignal);
            xxS = min(surround.generatorSignal) : max(surround.generatorSignal);
            rC = fitRes_indep.alphaC * ...
                normcdf(fitRes_indep.betaC * xxC + fitRes_indep.gammaC,0,1);
            rS = fitRes_indep.alphaS * ...
                normcdf(fitRes_indep.betaS.* xxS + fitRes_indep.gammaS,0,1);
            addLineToAxis(xxC,rC,...
                ['center',num2str(pp)],fig3,figColors(1,:),'-','none')
            addLineToAxis(xxS,rS,...
                ['surround',num2str(pp)],fig3,figColors(4,:),'-','none')
        end

        % 3) FIT SHARED NONLINEARITY MODEL. 5 FREE PARAMETERS
        % params is [a, alpha, beta, gamma, epsilon]
        params0=[2, 2*max(responseMean(:)), 0.4, 0, 0]';
        fitRes_shared = fitCSModel_SharedNL(centerGS,surroundGS,responseMean,params0);
            %2D fit surface:
            fitSurface = CSModel_SharedNL(CC(:)',SS(:)',...
                fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
                fitRes_shared.gamma,fitRes_shared.epsilon);
            fitSurface = reshape(fitSurface,length(ss),length(cc));

            subplot(223); hold off;
            stem3(centerGS,surroundGS,responseMean); hold on;
            surf(cc,ss,fitSurface)
            xlabel('Center'); ylabel('Surroud')
            title(['Shared, ',num2str(fitRes_shared.rSquared)])
            
        if currentNode.custom.get('isExample')
            xxCS = min(center.generatorSignal) + min(surround.generatorSignal) :...
                max(center.generatorSignal) + max(surround.generatorSignal);
            rCS = fitRes_shared.epsilon + fitRes_shared.alpha * ...
                normcdf(fitRes_shared.beta * xxCS + fitRes_shared.gamma,0,1);
            addLineToAxis(xxCS,rCS,...
                ['sharedNL',num2str(pp)],fig4,'k','-','none')
            
% % % % % % % % MODEL-FREE: SLICES THRU RESPONSE MATRIX % % % % % % % % % % % % % % %

            % Slices thru 2d nonlinearity
            colors = pmkmp(length(centerGS));
            for ii = 1:length(centerGS)
                addLineToAxis(surroundGS,responseMean(:,ii),...
                    ['centerSlice',num2str(ii)],fig6,colors(ii,:),'-','none')
                addLineToAxis(centerGS,responseMean(ii,:),...
                    ['surroundSlice',num2str(ii)],fig7,colors(ii,:),'-','none')
            end
            % C +/- Surround and vice-versa slices
            colors = pmkmp(length(centerGS));
            for ii = 1:length(centerGS)
                addLineToAxis(centerGS(ii) + surroundGS,responseMean(:,ii),...
                    ['centerRef',num2str(ii)],fig8,colors(ii,:),'-','none')
                
                addLineToAxis(surroundGS(ii) + centerGS,responseMean(ii,:),...
                    ['surroundRef',num2str(ii)],fig9,colors(ii,:),'-','none')
            end
            
            figure(17); clf; set(gcf, 'WindowStyle', 'docked')
            h = surf(centerGS,surroundGS,responseMean);
            h.FaceAlpha = 0.4;
            h.EdgeAlpha = 0.5;
            xlabel('Center activation','FontSize',14);
            ylabel('Surround activation','FontSize',14);
            zlabel('Response (nS)','FontSize',14)
            colormap(gray); caxis([-2, 15])
            hold on; plot3(centerGS,...
                surroundGS(1).*ones(size(centerGS)),...
                responseMean(1,:),'Color',colors(1,:),...
                'LineWidth',5,'Marker','none')
            hold on; plot3(centerGS,...
                surroundGS(end).*ones(size(centerGS)),...
                responseMean(end,:),'Color',colors(end,:),...
                'LineWidth',5,'Marker','none')
            view(-21,16)
            
            % Devation from linear summation
            linearSumResponse = centerMean + surroundMean;
            linearDeviationMatrix = responseMean - linearSumResponse;
            figure(18); clf;
            subplot(121);
            hs1 = surf(centerGS,surroundGS,responseMean); 
            colormap(hot); caxis([-2, 15]); freezeColors;
            xlabel('Center activation','FontSize',14);
            ylabel('Surround activation','FontSize',14);
            zlabel('Response (nS)','FontSize',14)
            hold on;
            hs2 = surf(centerGS,surroundGS,linearSumResponse);
            colormap(gray); caxis([-2, 15]); freezeColors;
            hs2.FaceAlpha = 0.4;
            hs2.EdgeAlpha = 0.5;
            view(-21,16)
            subplot(122);
            hs3 = surf(centerGS,surroundGS,linearDeviationMatrix);
            xlabel('Center activation','FontSize',14);
            ylabel('Surround activation','FontSize',14);
            zlabel('Measured - Linear Sum (nS)','FontSize',14)
            view(-21,16)
            
 
            %natural image luminances in 2D space
            load('~/Documents/MATLAB/turner-package/resources/SaccadeLuminanceTrajectoryStimuli_20160919.mat')
            numberOfBins_em = 100^2;
            centerGenSignal = [];
            surroundGenSignal = [];
            allCStim = [];
            allSStim = [];
            for ss = 1:length(luminanceData)
                cStim = resample(luminanceData(ss).centerTrajectory,center.sampleRate,200);
                cStim = (cStim) ./ luminanceData(ss).ImageMax; %stim as presented

                %convert to contrast (relative to mean) for filter convolution
                imMean = (luminanceData(ss).ImageMean  ./ luminanceData(ss).ImageMax);
                cStim = (cStim - imMean) / imMean;
                allCStim = cat(2,allCStim,cStim);

                linearPrediction = conv(cStim,center.LinearFilter);
                linearPrediction = linearPrediction(1:length(cStim));
                centerGenSignal = cat(2,centerGenSignal,linearPrediction);

                sStim = resample(luminanceData(ss).surroundTrajectory,surround.sampleRate,200);

                sStim = (sStim) ./ luminanceData(ss).ImageMax; %stim as presented

                %convert to contrast (relative to mean) for filter convolution
                imMean = (luminanceData(ss).ImageMean  ./ luminanceData(ss).ImageMax);
                sStim = (sStim - imMean) / imMean;
                allSStim = cat(2,allSStim,sStim);

                linearPrediction = conv(sStim,surround.LinearFilter);
                linearPrediction = linearPrediction(1:length(sStim));
                surroundGenSignal = cat(2,surroundGenSignal,linearPrediction);
            end

            figure; clf; fig19=gca; %2D stimulus space with eye movements
            set(gcf, 'WindowStyle', 'docked')
            set(fig19,'XScale','linear','YScale','linear')
            set(0, 'DefaultAxesFontSize', 12)
            set(get(fig19,'XLabel'),'String','Center gen. signal')
            set(get(fig19,'YLabel'),'String','Surround gen. signal')
            set(get(fig19,'YLabel'),'String','Probability')
            [N,Xedges,Yedges] = histcounts2(centerGenSignal,surroundGenSignal,sqrt(numberOfBins_em),...
                'Normalization','probability');
            Ccenters = Xedges(1:end-1) + diff(Xedges);
            Scenters = Yedges(1:end-1) + diff(Yedges);
            surf(Ccenters,Scenters,log10(N))
%             histogram2(centerGenSignal,surroundGenSignal,sqrt(numberOfBins_em),...
%                 'Normalization','probability','ShowEmptyBins','on','FaceColor','flat');
            colormap(hot)
            xlabel('Center'); ylabel('Surround'); zlabel('Probability')
            colorbar
            
        end

% % % % % % % % PREDICTIONS % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        ss_total = sum((testingData.csMeasured-mean(testingData.csMeasured)).^2);

        predictedResponse_joint = JointNLin_mvcn(testingData.centerGS, testingData.surroundGS,...
            fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);
        ss_resid_joint = sum((predictedResponse_joint'-testingData.csMeasured).^2);
        rSquaredValues.joint(pp) = 1-ss_resid_joint/ss_total;

        predictedResponse_indep = CSModel_IndependentNL(testingData.centerGS, testingData.surroundGS,...
            fitRes_indep.alphaC,fitRes_indep.betaC,...
            fitRes_indep.gammaC,fitRes_indep.alphaS,...
            fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);
        ss_resid_indep = sum((predictedResponse_indep-testingData.csMeasured).^2);
        rSquaredValues.indep(pp) = 1-ss_resid_indep/ss_total;

        predictedResponse_shared = CSModel_SharedNL(testingData.centerGS, testingData.surroundGS,...
            fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
            fitRes_shared.gamma,fitRes_shared.epsilon);
        ss_resid_shared = sum((predictedResponse_shared-testingData.csMeasured).^2);
        rSquaredValues.shared(pp) = 1-ss_resid_shared/ss_total;

    end
    
    addLineToAxis(rSquaredValues.indep(ONcellInds),rSquaredValues.shared(ONcellInds),...
        'ONr2',fig5,'b','none','o')
    addLineToAxis(rSquaredValues.indep(OFFcellInds),rSquaredValues.shared(OFFcellInds),...
        'OFFr2',fig5,'r','none','o')
    addLineToAxis([0 1],[0 1],...
        'unity',fig5,'k','--','none')
    disp(rSquaredValues);
    
    relativeImprovement = rSquaredValues.shared ./ rSquaredValues.indep;
    surroundWts = trapz(abs(filters.surround),2);
    centerWts = trapz(abs(filters.center),2);
    relativeSurroundWeight = surroundWts ./ centerWts;
    
    addLineToAxis(relativeSurroundWeight(ONcellInds),relativeImprovement(ONcellInds),...
        'ONimprovement',fig16,'b','none','o')
    addLineToAxis(relativeSurroundWeight(OFFcellInds),relativeImprovement(OFFcellInds),...
        'OFFimprovement',fig16,'r','none','o')
    addLineToAxis([0 1.1*max(relativeSurroundWeight)],[1 1],...
        'oneLine',fig16,'k','--','none')
    
    [rho, pval] = corr(relativeSurroundWeight,relativeImprovement');
    disp([rho, pval])
    
    recID = getRecordingTypeFromEpochList(currentNode.epochList);
    if (exportFigs)
        figID = ['CSLNfilters_',recID];
        makeAxisStruct(fig1,figID ,'RFSurroundFigs')

        figID = ['CSLNnls_',recID];
        makeAxisStruct(fig2,figID ,'RFSurroundFigs')

        figID = ['CSLNindNL_',recID];
        makeAxisStruct(fig3,figID ,'RFSurroundFigs')

        figID = ['CSLNsharedNL_',recID];
        makeAxisStruct(fig4,figID ,'RFSurroundFigs')

        figID = ['CSLNpopR2_',recID];
        makeAxisStruct(fig5,figID ,'RFSurroundFigs')

        figID = ['CSLNslice_C_',recID];
        makeAxisStruct(fig6,figID ,'RFSurroundFigs')

        figID = ['CSLNslice_S_',recID];
        makeAxisStruct(fig7,figID ,'RFSurroundFigs')

    %     figID = ['CSLNcs_',cellInfo.cellType,'_',recID];
    %     makeAxisStruct(fig8,figID ,'RFSurroundFigs')
    % 
    %     figID = ['CSLNsc_',cellInfo.cellType,'_',recID];
    %     makeAxisStruct(fig9,figID ,'RFSurroundFigs')

        figID = ['Cstim_',recID];
        makeAxisStruct(fig11,figID ,'RFSurroundFigs')

        figID = ['Sstim_',recID];
        makeAxisStruct(fig12,figID ,'RFSurroundFigs')

        figID = ['Cresp_',recID];
        makeAxisStruct(fig13,figID ,'RFSurroundFigs')

        figID = ['Sresp_',recID];
        makeAxisStruct(fig14,figID ,'RFSurroundFigs')

        figID = ['CSresp_',recID];
        makeAxisStruct(fig15,figID ,'RFSurroundFigs')
        
        figID = ['CSLNpopImprove_',recID];
        makeAxisStruct(fig16,figID ,'RFSurroundFigs')
    end
end