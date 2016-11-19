function doCSLNAnalysis(node,varargin)
    ip = inputParser;
    ip.addRequired('node',@(x)isa(x,'edu.washington.rieke.jauimodel.AuiEpochTree'));
    addParameter(ip,'bins2D',6^2,@isnumeric);
    addParameter(ip,'bins1D',20,@isnumeric); 
    
    ip.parse(node,varargin{:});
    node = ip.Results.node;
    bins2D = ip.Results.bins2D;
    bins1D = ip.Results.bins1D;
    
    figColors = pmkmp(8);

    figure; clf;
    fig1=gca; %Linear filters
    set(fig1,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig1,'XLabel'),'String','Time (s)')
    set(get(fig1,'YLabel'),'String','')
    
    figure; clf;
    fig2=gca; %Nonlinearities
    set(fig2,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig2,'XLabel'),'String','Linear prediction (pA)')
    set(get(fig2,'YLabel'),'String','Measured (pA)')
    
    figure; clf;
    fig3=gca; %Independent nonlinearity model
    set(fig3,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig3,'XLabel'),'String','Linear prediction (pA)')
    set(get(fig3,'YLabel'),'String','Measured (pA)')
    
    figure; clf;
    fig4=gca; %Shared nonlinearity
    set(fig4,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig4,'XLabel'),'String','Linear prediction (pA)')
    set(get(fig4,'YLabel'),'String','Measured (pA)')
    
    figure; clf;
    fig5=gca; %population R-squared results
    set(fig5,'XScale','linear','YScale','linear')
    set(0, 'DefaultAxesFontSize', 12)
    set(get(fig5,'XLabel'),'String','R^2 independent nonlinearities')
    set(get(fig5,'YLabel'),'String','R^2 shared nonlinearity')


    populationNodes = {};
    ct = 0;
    for nn = 1:node.descendentsDepthFirst.length
        if strcmp(node.descendentsDepthFirst(nn).splitKey,...
                'protocolSettings(currentStimulus)') && node.descendentsDepthFirst(nn).custom.get('isSelected')
            ct = ct + 1;
            populationNodes(ct) = node.descendentsDepthFirst(nn); %#ok<AGROW>
        end
    end
    
    rSquaredValues.joint = [];
    rSquaredValues.indep = [];
    rSquaredValues.shared = [];
    for pp = 1:length(populationNodes)
        currentNode = populationNodes{pp};
        cellInfo = getCellInfoFromEpochList(currentNode.epochList);
        recType = getRecordingTypeFromEpochList(currentNode.epochList);
        
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

        if currentNode.custom.get('isExample')
                addLineToAxis(center.filterTimeVector,center.LinearFilter,...
                    ['center',num2str(pp)],fig1,figColors(1,:),'-','none')
                addLineToAxis(surround.filterTimeVector,surround.LinearFilter,...
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
        end

% % % % % % % % MODEL FITTING % % % % % % % % % % % % % % % % % % % % % % % % % 
        epochLen = length(centerSurround.measuredResponse) / centerSurround.n;
        epochToHold = 5; %epoch number to hold out
        testDataInds = ((epochToHold-1)*epochLen + 1):(epochToHold*epochLen);
        fitDataInds = setdiff(1:length(centerSurround.measuredResponse),testDataInds);

        %bin up and shape training data:
        [~,centerGS,~,centerBinID] = ...
            histcounts_equallyPopulatedBins(center.generatorSignal(fitDataInds),sqrt(bins2D));
        [~,surroundGS,~,surroundBinID] = ...
            histcounts_equallyPopulatedBins(surround.generatorSignal(fitDataInds),sqrt(bins2D));
        responseMean = zeros(sqrt(bins2D));
        for xx = 1:sqrt(bins2D)
            for yy = 1:sqrt(bins2D)
                jointInds = intersect(find(centerBinID == xx),find(surroundBinID == yy));
                tempResp = centerSurround.measuredResponse(fitDataInds);
                responseMean(yy,xx) = mean(tempResp(jointInds)); %responseMean is [surround, center]
            end
        end
        %1) FIT JOINT NONLINEARITY MODEL. 7 FREE PARAMETERS
        % params = [alpha mu1 mu2 std1 std2 corr12 epsilon]
        params0 = [max(responseMean(:)), 400, 0, 50, 50, 0, 0];
        fitRes_joint = fitNLinearity_2D(centerGS,surroundGS,responseMean,params0);
            %2D fit surface:
            cc = linspace(min(centerGS),max(centerGS),20);
            ss = linspace(min(surroundGS),max(surroundGS),20);
            [CC,SS] = meshgrid(cc',ss');
            fitSurface = JointNLin_mvcn(CC(:)',SS(:)',fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);
            fitSurface = reshape(fitSurface,length(ss),length(cc)); %#ok<NASGU>
        
        % 2) FIT INDEPENDENT NONLINEARITY MODEL. 7 FREE PARAMETERS
        %params is [alphaC, betaC, gammaC,...
        %           alphaS, betaS, gammaS, epsilon]
        params0 = [max(responseMean(:)), center.nonlinearity.fitParams.beta, center.nonlinearity.fitParams.gamma,...
            max(responseMean(:)), surround.nonlinearity.fitParams.beta, surround.nonlinearity.fitParams.gamma,...
            min([surround.nonlinearity.fitParams.epsilon, center.nonlinearity.fitParams.epsilon])];
        fitRes_indep = fitCSModel_IndependentNL(centerGS,surroundGS,responseMean,params0);
            %2D fit surface:
            fitSurface = CSModel_IndependentNL(CC(:)',SS(:)',...
                fitRes_indep.alphaC,fitRes_indep.betaC,...
                fitRes_indep.gammaC,fitRes_indep.alphaS,...
                fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);

            fitSurface = reshape(fitSurface,length(ss),length(cc)); %#ok<NASGU>
            
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
        params0=[2, max(responseMean(:)), 0.01, 0, 0]';
        fitRes_shared = fitCSModel_SharedNL(centerGS,surroundGS,responseMean,params0);
            %2D fit surface:
            fitSurface = CSModel_SharedNL(CC(:)',SS(:)',...
                fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
                fitRes_shared.gamma,fitRes_shared.epsilon);
            fitSurface = reshape(fitSurface,length(ss),length(cc)); %#ok<NASGU>
            
        if currentNode.custom.get('isExample')
            xxCS = min(center.generatorSignal) + min(surround.generatorSignal) :...
                max(center.generatorSignal) + max(surround.generatorSignal);
            rCS = fitRes_shared.epsilon + fitRes_shared.alpha * ...
                normcdf(fitRes_shared.beta * xxCS + fitRes_shared.gamma,0,1);
            addLineToAxis(xxCS,rCS,...
                ['sharedNL',num2str(pp)],fig4,'k','-','none')
            
% % % % % % % % MODEL-FREE: SLICES THRU RESPONSE MATRIX % % % % % % % % % % % % % % %
            figure; clf;
            fig6=gca; %slice of center, modulate surround
            set(fig6,'XScale','linear','YScale','linear')
            set(0, 'DefaultAxesFontSize', 12)
            set(get(fig6,'XLabel'),'String','Surround activation')
            set(get(fig6,'YLabel'),'String','Response (pA)')

            figure; clf;
            fig7=gca; %slice of surround, modulate center
            set(fig7,'XScale','linear','YScale','linear')
            set(0, 'DefaultAxesFontSize', 12)
            set(get(fig7,'XLabel'),'String','Center activation')
            set(get(fig7,'YLabel'),'String','Response (pA)')

            colors = pmkmp(length(centerGS));
            for ii = 1:length(centerGS)
                addLineToAxis(surroundGS,responseMean(:,ii),...
                    ['centerSlice',num2str(ii)],fig6,colors(ii,:),'-','none')
                addLineToAxis(centerGS,responseMean(ii,:),...
                    ['surroundSlice',num2str(ii)],fig7,colors(ii,:),'-','none')
            end
            
            figure; clf;
            fig8=gca; %center +/- surround overlay
            set(fig8,'XScale','linear','YScale','linear')
            set(0, 'DefaultAxesFontSize', 12)
            set(get(fig8,'XLabel'),'String','Center +/- surround')
            set(get(fig8,'YLabel'),'String','Response (pA)')
            
            figure; clf;
            fig9=gca; %surround +/- center overlay
            set(fig9,'XScale','linear','YScale','linear')
            set(0, 'DefaultAxesFontSize', 12)
            set(get(fig9,'XLabel'),'String','Surround +/- center')
            set(get(fig9,'YLabel'),'String','Response (pA)')
            
            colors = pmkmp(length(centerGS));
            for ii = 1:length(centerGS)
                addLineToAxis(centerGS(ii) + surroundGS,responseMean(:,ii),...
                    ['centerRef',num2str(ii)],fig8,colors(ii,:),'-','none')
                
                addLineToAxis(surroundGS(ii) + centerGS,responseMean(ii,:),...
                    ['surroundRef',num2str(ii)],fig9,colors(ii,:),'-','none')
            end
            
            figure; clf;
            fig10=gca; %2D mesh response matrix
            set(fig10,'XScale','linear','YScale','linear')
            set(0, 'DefaultAxesFontSize', 12)
            set(get(fig10,'XLabel'),'String','Center')
            set(get(fig10,'XLabel'),'String','Surround')
            set(get(fig10,'YLabel'),'String','Response (pA)')
            surfc(centerGS,surroundGS,responseMean);
            xlabel('Center'); ylabel('Surround'); zlabel('Response (pA)')
            
        end

% % % % % % % % PREDICTIONS % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        centerActivation = center.generatorSignal(testDataInds);
        surroundActivation = surround.generatorSignal(testDataInds);
        measuredResponse = centerSurround.measuredResponse(testDataInds);
        ss_total = sum((measuredResponse-mean(measuredResponse)).^2);

        predictedResponse_joint = JointNLin_mvcn(centerActivation, surroundActivation,...
            fitRes_joint.alpha,fitRes_joint.mu,fitRes_joint.sigma,fitRes_joint.epsilon);
        ss_resid_joint = sum((predictedResponse_joint'-measuredResponse).^2);
        rSquaredValues.joint(pp) = 1-ss_resid_joint/ss_total;

        predictedResponse_indep = CSModel_IndependentNL(centerActivation, surroundActivation,...
            fitRes_indep.alphaC,fitRes_indep.betaC,...
            fitRes_indep.gammaC,fitRes_indep.alphaS,...
            fitRes_indep.betaS,fitRes_indep.gammaS,fitRes_indep.epsilon);
        ss_resid_indep = sum((predictedResponse_indep-measuredResponse).^2);
        rSquaredValues.indep(pp) = 1-ss_resid_indep/ss_total;

        predictedResponse_shared = CSModel_SharedNL(centerActivation, surroundActivation,...
            fitRes_shared.a,fitRes_shared.alpha,fitRes_shared.beta,...
            fitRes_shared.gamma,fitRes_shared.epsilon);
        ss_resid_shared = sum((predictedResponse_shared-measuredResponse).^2);
        rSquaredValues.shared(pp) = 1-ss_resid_shared/ss_total;
    end
    
    addLineToAxis(rSquaredValues.indep,rSquaredValues.shared,...
        'r2',fig5,'k','none','o')
    addLineToAxis([0 1],[0 1],...
        'unity',fig5,'k','--','none')
    
    
    figID = ['CSLNfilters_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig1,figID ,'RFSurroundFigs')
    
    figID = ['CSLNnls_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig2,figID ,'RFSurroundFigs')
    
    figID = ['CSLNindNL_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig3,figID ,'RFSurroundFigs')
    
    figID = ['CSLNsharedNL_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig4,figID ,'RFSurroundFigs')
    
    figID = ['CSLNpopR2_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig5,figID ,'RFSurroundFigs')
    
    figID = ['CSLNslice_C_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig6,figID ,'RFSurroundFigs')

    figID = ['CSLNslice_S_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig7,figID ,'RFSurroundFigs')
    
    figID = ['CSLNcs_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig8,figID ,'RFSurroundFigs')

    figID = ['CSLNsc_',cellInfo.cellType,'_',recType];
    makeAxisStruct(fig9,figID ,'RFSurroundFigs')

end