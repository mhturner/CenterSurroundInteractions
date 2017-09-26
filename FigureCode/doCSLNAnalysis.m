function doCSLNAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'bins2D',6^2,@isnumeric);
    addParameter(ip,'bins1D',20,@isnumeric); 
    addParameter(ip,'exportFigs',true,@islogical);
    addParameter(ip,'convertToConductance',true,@islogical);
    addParameter(ip,'fitWithEquallyPopulatedBins',true,@islogical);
    
    figDir = '~/Dropbox/RiekeLab/Analysis/MATLAB/RFSurround/resources/TempFigs/'; %for saved eps figs

    ip.parse(node,varargin{:});
    node = ip.Results.node;
    bins2D = ip.Results.bins2D;
    bins1D = ip.Results.bins1D;
    exportFigs = ip.Results.exportFigs;
    convertToConductance = ip.Results.convertToConductance;
    fitWithEquallyPopulatedBins = ip.Results.fitWithEquallyPopulatedBins;
    
    figColors = pmkmp(8);
    egEpochToPull = 3; %for example traces
    
    figure; clf; fig1=gca; initFig(fig1,'Time (s)','') % Linear filters
    figure; clf; fig2=gca; initFig(fig2,'Linear prediction (nS)','Measured (nS)') % Nonlinearities
    figure; clf; fig3=gca; initFig(fig3,'Linear prediction (nS)','Measured (nS)') % Independent nonlinearities
    figure; clf; fig4=gca; initFig(fig4,'Linear prediction (nS)','Measured (nS)') % Shared nonlinearity
    
    figure; clf; fig5=gca; initFig(fig5,'R^2 independent nonlinearities','R^2 shared nonlinearity') % pop R2. Shared vs indep
    figure; clf; fig10=gca; initFig(fig10,'R^2 Three nonlinearities','R^2 shared nonlinearity') % pop R2. Shared vs indep
    
    figure; clf; fig6=gca; initFig(fig6,'Surround activation','Response (nS)') %slice of center, modulate surround
    figure; clf; fig7=gca; initFig(fig7,'Center activation','Response (nS)') %slice of surround, modulate center
   
    % Individual traces for example cell:
    figure; clf; fig8=gca; initFig(fig8,'Time (s)','Response (nS)') %Eg trace of R(C+S) and [R(C) + R(S)] overlay
    figure; clf; fig11=gca; initFig(fig11,'Time (s)','Contrast') %center stim
    figure; clf; fig12=gca; initFig(fig12,'Time (s)','Contrast') %surround stim
    figure; clf; fig13=gca; initFig(fig13,'Time (s)','Response (nS)') %center response
    figure; clf; fig14=gca; initFig(fig14,'Time (s)','Response (nS)') %surround response
    figure; clf; fig15=gca; initFig(fig15,'Time (s)','Response (nS)') %center-surround response

    figure; clf; fig16=gca; initFig(fig16,'Relative surround weight','Improvement over independent model') %improvement in R2 vs surround weight
    
    figure; clf; fig17=gca; initFig(fig17,'Linear prediction (nS)','Output (a.u.)') %ThreeNL model c & s nonlinearities
    figure; clf; fig18=gca; initFig(fig18,'Combined C + S input (a.u.)','Measured (nS)') %ThreeNL model output nonlinearity

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
    rSquaredValues.threeNL = [];
    ONcellInds = [];
    OFFcellInds = [];
    for pp = 1:length(populationNodes)
        recNode = populationNodes{pp};
        currentNode = [];
        repeatedNode = [];
        for cc = 1:recNode.children.length
           if recNode.children(cc).splitValue == 1
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
            %Linear filters
            addLineToAxis([center.filterTimeVector],[center.LinearFilter],...
                ['center',num2str(pp)],fig1,figColors(1,:),'-','none')
            addLineToAxis([surround.filterTimeVector],[surround.LinearFilter],...
                ['surround',num2str(pp)],fig1,figColors(4,:),'-','none')
            %Independently fit LN model nonlinearities (not CS model):
            %   ...Center
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
            %   ...Surround
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

            % Example Traces:
            %   ...Center stim and response
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
            %   ...Surround stim and response
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
            %   ...C + S response
            csEL = sortEpochList_time(currentNode.childBySplitValue('Center-Surround').epochList);
            csExampleEpoch = csEL.elements(egEpochToPull);
            cscSeed = csExampleEpoch.protocolSettings('centerNoiseSeed');
            cssSeed = csExampleEpoch.protocolSettings('surroundNoiseSeed');
            csEpochRes = getNoiseStimulusAndResponse(csExampleEpoch,recType,'seedName','centerNoiseSeed');
            addLineToAxis(csEpochRes.wholeTrace.timeVector(1:length(csEpochRes.wholeTrace.stimulus)),...
                csEpochRes.wholeTrace.response(1:length(csEpochRes.wholeTrace.stimulus)),...
                ['resp',num2str(pp)],fig15,'k','-','none')
            
            % Example Traces:
            %   C + S and response vs linear sum.
            %   i.e. R(C+S) overlay [R(C) + R(S)]
            tempCStruct = getMeanResponseTrace(repeatedNode.childBySplitValue('Center').epochList,recType);
            tempSStruct = getMeanResponseTrace(repeatedNode.childBySplitValue('Surround').epochList,recType);
            tempLinSum = tempCStruct.mean + tempSStruct.mean;
            
            tempStruct = getMeanResponseTrace(repeatedNode.childBySplitValue('Center-Surround').epochList,recType);

            addLineToAxis(tempStruct.timeVector,tempStruct.mean,'measuredResp',fig8,'k','-','none')
            addLineToAxis(tempStruct.timeVector,tempLinSum,'linSum',fig8,[0.8 0.8 0.8],'-','none')
 
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
            temp = [];
            for testingEpoch = 1:(centerSurround.n)
                epochToHold = testingEpoch; %epoch number to hold out
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
                
                predRes = getPredictions(trainingData,testingData,...
                    bins2D,fitWithEquallyPopulatedBins,...
                    center,surround,...
                    currentNode.custom.get('isExample'));
                temp(testingEpoch,1) = predRes.joint;
                temp(testingEpoch,2) = predRes.indep;
                temp(testingEpoch,3) = predRes.shared;
                temp(testingEpoch,4) = predRes.threeNL;
            end
            
            %average over all testing epochs:
            rSquaredValues.joint(pp) = mean(temp(:,1));
            rSquaredValues.indep(pp) = mean(temp(:,2));
            rSquaredValues.shared(pp) = mean(temp(:,3));
            rSquaredValues.threeNL(pp) = mean(temp(:,4));
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
            
            predRes = getPredictions(trainingData,testingData,...
                bins2D,fitWithEquallyPopulatedBins,...
                center,surround,...
                currentNode.custom.get('isExample'));
            rSquaredValues.joint(pp) = predRes.joint;
            rSquaredValues.indep(pp) = predRes.indep;
            rSquaredValues.shared(pp) = predRes.shared;
            rSquaredValues.threeNL(pp) = predRes.threeNL;
        end
        
% % % % % % % % MODEL-FREE: SLICES THRU RESPONSE MATRIX % % % % % % % % % % % % % % %
        if currentNode.custom.get('isExample')
            binRes = binData(trainingData,bins2D,fitWithEquallyPopulatedBins);
            centerGS = binRes.centerGS;
            surroundGS = binRes.surroundGS;
            responseMean = binRes.responseMean;
            
            % Slices thru 2d nonlinearity
            colors = pmkmp(length(centerGS));
            for ii = 1:length(centerGS)
                addLineToAxis(surroundGS,responseMean(:,ii),...
                    ['centerSlice',num2str(ii)],fig6,colors(ii,:),'-','none')
                addLineToAxis(centerGS,responseMean(ii,:),...
                    ['surroundSlice',num2str(ii)],fig7,colors(ii,:),'-','none')
            end

            %3D response surface with highlighted slices
            % write out to csv for igor..
            [tempCC, tempSS] = meshgrid(centerGS,surroundGS);
            csvwrite([figDir, 'Layer1.csv'],tempCC);
            csvwrite([figDir, 'Layer2.csv'],tempSS);
            csvwrite([figDir, 'Layer3.csv'],responseMean);

            save([figDir, 'Layer1.mat'],'tempCC');
            save([figDir, 'Layer2.mat'],'tempSS');
            save([figDir, 'Layer3.mat'],'responseMean');
                 
            %make in matlab..
            fh = figure(22); clf;
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
            set(fh,'Position',[1397         676         501         419])
            drawnow;
            figID = 'respSurf_highlightedSlices';
            print(fh,[figDir,figID],'-depsc')
        end %example plotting of model-free additivity stuff

    end
    
    %Population data, indep vs. shared:
    %   On...
    addLineToAxis(rSquaredValues.indep(ONcellInds),rSquaredValues.shared(ONcellInds),...
        'ONr2',fig5,'b','none','o')
    meanON_x = mean(rSquaredValues.indep(ONcellInds));
    meanON_y = mean(rSquaredValues.shared(ONcellInds));
    errON_x = std(rSquaredValues.indep(ONcellInds)) ./ sqrt(length(ONcellInds));
    errON_y = std(rSquaredValues.shared(ONcellInds)) ./ sqrt(length(ONcellInds));
    addLineToAxis(meanON_x,meanON_y,'ONmean',fig5,'b','none','.')
    addLineToAxis(meanON_x + [errON_x, -errON_x],[meanON_y, meanON_y],'ONerr_x',fig5,'b','-','none')
    addLineToAxis([meanON_x, meanON_x],meanON_y + [errON_y, -errON_y],'ONerr_y',fig5,'b','-','none')
    [~,p] = ttest(rSquaredValues.indep(ONcellInds), rSquaredValues.shared(ONcellInds));
    disp('Shared vs. Indep. On:')
    disp(['p = ',num2str(p)])
    
    %   Off...
    addLineToAxis(rSquaredValues.indep(OFFcellInds),rSquaredValues.shared(OFFcellInds),...
        'OFFr2',fig5,'r','none','o')
    meanOFF_x = mean(rSquaredValues.indep(OFFcellInds));
    meanOFF_y = mean(rSquaredValues.shared(OFFcellInds));
    errOFF_x = std(rSquaredValues.indep(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    errOFF_y = std(rSquaredValues.shared(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    addLineToAxis(meanOFF_x,meanOFF_y,'OFFmean',fig5,'r','none','.')
    addLineToAxis(meanOFF_x + [errOFF_x, -errOFF_x],[meanOFF_y, meanOFF_y],'OFFerr_x',fig5,'r','-','none')
    addLineToAxis([meanOFF_x, meanOFF_x],meanOFF_y + [errOFF_y, -errOFF_y],'OFFerr_y',fig5,'r','-','none')
    [~, p] = ttest(rSquaredValues.indep(OFFcellInds), rSquaredValues.shared(OFFcellInds));
    disp('Shared vs. Indep. Off:')
    disp(['p = ',num2str(p)])
    
    addLineToAxis([0 1],[0 1],...
        'unity',fig5,'k','--','none')

    %Population data, threeNL (stacked) vs. shared:
    %   On...
    addLineToAxis(rSquaredValues.threeNL(ONcellInds),rSquaredValues.shared(ONcellInds),...
        'ONr2',fig10,'b','none','o')
    meanON_x = mean(rSquaredValues.threeNL(ONcellInds));
    meanON_y = mean(rSquaredValues.shared(ONcellInds));
    errON_x = std(rSquaredValues.threeNL(ONcellInds)) ./ sqrt(length(ONcellInds));
    errON_y = std(rSquaredValues.shared(ONcellInds)) ./ sqrt(length(ONcellInds));
    addLineToAxis(meanON_x,meanON_y,'ONmean',fig10,'b','none','.')
    addLineToAxis(meanON_x + [errON_x, -errON_x],[meanON_y, meanON_y],'ONerr_x',fig10,'b','-','none')
    addLineToAxis([meanON_x, meanON_x],meanON_y + [errON_y, -errON_y],'ONerr_y',fig10,'b','-','none')
    [~, p] = ttest(rSquaredValues.threeNL(ONcellInds), rSquaredValues.shared(ONcellInds));
    disp('Shared vs. ThreeNL On:')
    disp(['p = ',num2str(p)])
    
    addLineToAxis(rSquaredValues.threeNL(OFFcellInds),rSquaredValues.shared(OFFcellInds),...
        'OFFr2',fig10,'r','none','o')
    meanOFF_x = mean(rSquaredValues.threeNL(OFFcellInds));
    meanOFF_y = mean(rSquaredValues.shared(OFFcellInds));
    errOFF_x = std(rSquaredValues.threeNL(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    errOFF_y = std(rSquaredValues.shared(OFFcellInds)) ./ sqrt(length(OFFcellInds));
    addLineToAxis(meanOFF_x,meanOFF_y,'OFFmean',fig10,'r','none','.')
    addLineToAxis(meanOFF_x + [errOFF_x, -errOFF_x],[meanOFF_y, meanOFF_y],'OFFerr_x',fig10,'r','-','none')
    addLineToAxis([meanOFF_x, meanOFF_x],meanOFF_y + [errOFF_y, -errOFF_y],'OFFerr_y',fig10,'r','-','none')
    [~, p] = ttest(rSquaredValues.threeNL(OFFcellInds), rSquaredValues.shared(OFFcellInds));
    disp('Shared vs. ThreeNL Off:')
    disp(['p = ',num2str(p)])
    
    addLineToAxis([0 1],[0 1],...
        'unity',fig10,'k','--','none')
    
    %improvement vs. surround weight
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
    
    [rho, pval] = corr(relativeSurroundWeight(ONcellInds),relativeImprovement(ONcellInds)');
    disp('Surr. weight corr with improvement, ON:')
    disp([rho, pval])
    
    [rho, pval] = corr(relativeSurroundWeight(OFFcellInds),relativeImprovement(OFFcellInds)');
    disp('Surr. weight corr with improvement, OFF:')
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

        figID = ['CSLNlinSumTrace_','_',recID];
        makeAxisStruct(fig8,figID ,'RFSurroundFigs')
    
        figID = ['CSLNpopR2_3NL_',recID];
        makeAxisStruct(fig10,figID ,'RFSurroundFigs')

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
        
        figID = ['CSLN_3NLcs',recID];
        makeAxisStruct(fig17,figID ,'RFSurroundFigs')
        
        figID = ['CSLN_3NLshared',recID];
        makeAxisStruct(fig18,figID ,'RFSurroundFigs')
    end
end

function res = binData(data,bins2D,fitWithEquallyPopulatedBins)
 %bin up and shape training data:
    if (fitWithEquallyPopulatedBins)
        %equally populated bins...
        [~,res.centerGS,~,centerBinID] = ...
            histcounts_equallyPopulatedBins(data.centerGS,sqrt(bins2D));
        [~,res.surroundGS,~,surroundBinID] = ...
            histcounts_equallyPopulatedBins(data.surroundGS,sqrt(bins2D));
    else
        %evenly-spaced bins...
        [~,edges,centerBinID] = histcounts(data.centerGS,sqrt(bins2D));
        res.centerGS = edges(1:end-1) + diff(edges);
        [~,edges,surroundBinID] = histcounts(data.surroundGS,sqrt(bins2D));
        res.surroundGS = edges(1:end-1) + diff(edges);
    end

    res.centerMean = zeros(sqrt(bins2D));
    res.surroundMean = zeros(sqrt(bins2D));
    res.responseMean = zeros(sqrt(bins2D));
    res.responseErr = zeros(sqrt(bins2D));
    for xx = 1:sqrt(bins2D)
    for yy = 1:sqrt(bins2D)
        jointInds = intersect(find(centerBinID == xx),find(surroundBinID == yy));
        res.centerMean(yy,xx) = mean(data.cMeasured(jointInds));
        res.surroundMean(yy,xx) = mean(data.sMeasured(jointInds));
        res.responseMean(yy,xx) = mean(data.csMeasured(jointInds));
        res.responseErr(yy,xx) = std(data.csMeasured(jointInds));
    end
    end

end

function res = getPredictions(trainingData,testingData,bins2D,...
fitWithEquallyPopulatedBins,center,surround,isExample)
    figColors = pmkmp(8);
    binRes = binData(trainingData,bins2D,fitWithEquallyPopulatedBins);
    centerGS = binRes.centerGS;
    surroundGS = binRes.surroundGS;
    responseMean = binRes.responseMean;

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

    figure(21); clf; set(gcf, 'WindowStyle', 'docked')
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
%     if pp == 7
%     betaGuess = 0.1;
%     upGuess = 2*max(responseMean(:));
%     end
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

    figure(21);
    subplot(222); hold off;
    stem3(centerGS,surroundGS,responseMean); hold on;
    surf(cc,ss,fitSurface)
    xlabel('Center'); ylabel('Surroud')
    title(['Indep, ',num2str(fitRes_indep.rSquared)])

%     if (isExample)
%     xxC = min(center.generatorSignal) : max(center.generatorSignal);
%     xxS = min(surround.generatorSignal) : max(surround.generatorSignal);
%     rC = fitRes_indep.alphaC * ...
%         normcdf(fitRes_indep.betaC * xxC + fitRes_indep.gammaC,0,1);
%     rS = fitRes_indep.alphaS * ...
%         normcdf(fitRes_indep.betaS.* xxS + fitRes_indep.gammaS,0,1);
%     addLineToAxis(xxC,rC,...
%         ['center'],fig3,figColors(1,:),'-','none')
%     addLineToAxis(xxS,rS,...
%         ['surround'],fig3,figColors(4,:),'-','none')
%     end

    % 3) FIT SHARED NONLINEARITY MODEL. 5 FREE PARAMETERS
    % params is [a, alpha, beta, gamma, epsilon]
    params0=[2, 2*max(responseMean(:)), 0.4, 0, 0]';
    fitRes_shared = fitCSModel_SharedNL(centerGS,surroundGS,responseMean,params0);
    %2D fit surface:
    fitSurface = CSModel_SharedNL(CC(:)',SS(:)',...
        fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
        fitRes_shared.gamma,fitRes_shared.epsilon);
    fitSurface = reshape(fitSurface,length(ss),length(cc));

    figure(21);
    subplot(223); hold off;
    stem3(centerGS,surroundGS,responseMean); hold on;
    surf(cc,ss,fitSurface)
    xlabel('Center'); ylabel('Surroud')
    title(['Shared, ',num2str(fitRes_shared.rSquared)])

%     if (isExample)
%         xxCS = min(center.generatorSignal) + min(surround.generatorSignal) :...
%             max(center.generatorSignal) + max(surround.generatorSignal);
%         rCS = fitRes_shared.epsilon + fitRes_shared.alpha * ...
%             normcdf(fitRes_shared.beta * xxCS + fitRes_shared.gamma,0,1);
%         addLineToAxis(xxCS,rCS,...
%             ['sharedNL'],fig4,'k','-','none')
%     end


    % 4) FIT THREE NONLINEARITY MODEL. 10 FREE PARAMETERS
    %params is [alphaC, betaC, gammaC,...
    %           alphaS, betaS, gammaS, ...
    %           alphaShared, betaShared, gammaShared, epsilon]
    upGuess = 3*max(responseMean(:));
    betaGuess = center.nonlinearity.fitParams.beta;
    epsilonGuess  = min([surround.nonlinearity.fitParams.epsilon, center.nonlinearity.fitParams.epsilon]);

    params0 = [1, 0.1, 0,...
    1, 0.1, 0,...
    upGuess, betaGuess, center.nonlinearity.fitParams.gamma,...
    epsilonGuess];
    fitRes_ThreeNL = fitCSModel_ThreeNL(centerGS,surroundGS,responseMean,params0);
    %2D fit surface:
    fitSurface = CSModel_ThreeNL(CC(:)',SS(:)',...
        fitRes_ThreeNL.alphaC,fitRes_ThreeNL.betaC,fitRes_ThreeNL.gammaC,...
        fitRes_ThreeNL.alphaS,fitRes_ThreeNL.betaS,fitRes_ThreeNL.gammaS,...
        fitRes_ThreeNL.alphaShared,fitRes_ThreeNL.betaShared,fitRes_ThreeNL.gammaShared,...
        fitRes_ThreeNL.epsilon);

    fitSurface = reshape(fitSurface,length(ss),length(cc));

    figure(21);
    subplot(224); hold off;
    stem3(centerGS,surroundGS,responseMean); hold on;
    surf(cc,ss,fitSurface)
    xlabel('Center'); ylabel('Surroud')
    title(['ThreeNL, ',num2str(fitRes_ThreeNL.rSquared)])

%     if (isExample)
%         xxC = min(center.generatorSignal) : max(center.generatorSignal);
%         xxS = min(surround.generatorSignal) : max(surround.generatorSignal);
%         rC = fitRes_ThreeNL.alphaC * ...
%             normcdf(fitRes_ThreeNL.betaC * xxC + fitRes_ThreeNL.gammaC,0,1);
%         rS = fitRes_ThreeNL.alphaS * ...
%             normcdf(fitRes_ThreeNL.betaS.* xxS + fitRes_ThreeNL.gammaS,0,1);
%         xxCS = linspace(min(rC) + min(rS), max(rC) + max(rS),100);
%         rOut = fitRes_ThreeNL.alphaShared * ...
%             normcdf(fitRes_ThreeNL.betaShared.* xxCS + fitRes_ThreeNL.gammaShared,0,1);
%         addLineToAxis(xxC,rC,...
%             ['center'],fig17,figColors(1,:),'-','none')
%         addLineToAxis(xxS,rS,...
%             ['surround'],fig17,figColors(4,:),'-','none')
%         addLineToAxis(xxCS,rOut,...
%             ['output'],fig18,'k','-','none')
%         cZero = fitRes_ThreeNL.alphaC * ...
%             normcdf(fitRes_ThreeNL.betaC * 0 + fitRes_ThreeNL.gammaC,0,1);
%         sZero = fitRes_ThreeNL.alphaS * ...
%             normcdf(fitRes_ThreeNL.betaS.* 0 + fitRes_ThreeNL.gammaS,0,1);
%         baselineXlevel =  cZero + sZero;
%         addLineToAxis([baselineXlevel baselineXlevel],[0 max(rOut)],...
%             ['baselineX',num2str(pp)],fig18,'k','--','none')
% 
%     end
    % % % % % % % % PREDICTIONS % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    ss_total = sum((testingData.csMeasured-mean(testingData.csMeasured)).^2);

    predictedResponse_joint = JointNLin_mvcn(testingData.centerGS, testingData.surroundGS,...
        fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);
    ss_resid_joint = sum((predictedResponse_joint'-testingData.csMeasured).^2);
    res.joint = 1-ss_resid_joint/ss_total;

    predictedResponse_indep = CSModel_IndependentNL(testingData.centerGS, testingData.surroundGS,...
        fitRes_indep.alphaC,fitRes_indep.betaC,...
        fitRes_indep.gammaC,fitRes_indep.alphaS,...
        fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);
    ss_resid_indep = sum((predictedResponse_indep-testingData.csMeasured).^2);
    res.indep = 1-ss_resid_indep/ss_total;

    predictedResponse_shared = CSModel_SharedNL(testingData.centerGS, testingData.surroundGS,...
        fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
        fitRes_shared.gamma,fitRes_shared.epsilon);
    ss_resid_shared = sum((predictedResponse_shared-testingData.csMeasured).^2);
    res.shared = 1-ss_resid_shared/ss_total;

    predictedResponse_threeNL = CSModel_ThreeNL(testingData.centerGS, testingData.surroundGS,...
        fitRes_ThreeNL.alphaC,fitRes_ThreeNL.betaC,fitRes_ThreeNL.gammaC,...
        fitRes_ThreeNL.alphaS,fitRes_ThreeNL.betaS,fitRes_ThreeNL.gammaS,...
        fitRes_ThreeNL.alphaShared,fitRes_ThreeNL.betaShared,fitRes_ThreeNL.gammaShared,...
        fitRes_ThreeNL.epsilon);
    ss_resid_ThreeNL = sum((predictedResponse_threeNL-testingData.csMeasured).^2);
    res.threeNL = 1-ss_resid_ThreeNL/ss_total;


end