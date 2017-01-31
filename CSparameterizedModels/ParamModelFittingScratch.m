



timeVector = linspace(0,1,1e3);
linearFilter = ParameterizedLinearFilter(timeVector,'Phase',-10);

stimTime = 0.001:0.001:2;
stimulus = zeros(size(stimTime));
stimulus(100:300) = 1;
stimulus(600:800) = -1;

linearResponse = conv(stimulus,linearFilter);
linearResponse = linearResponse(1:length(stimTime));

figure(1); clf; subplot(3,2,1); plot(timeVector,linearFilter,'b')
subplot(3,2,3:4); plot(stimTime,stimulus);
subplot(3,2,5:6); plot(stimTime,linearResponse,'k');


%%

response = ParameterizedNonlinearity(linearResponse);

figure(1); clf; plot(timeVector,linearResponse,'b')
figure(1); hold on; plot(timeVector,response,'k')


%%

centerFilterParams.NumFilters = 10;
centerFilterParams.TauRise = 0.1;
centerFilterParams.TauDamp = 0.03;
centerFilterParams.Period = 0.5;
centerFilterParams.Phase = -10; % % %
centerFilterParams.ScaleFactor = 10;

surroundFilterParams.NumFilters = 10;
surroundFilterParams.TauRise = 0.1;
surroundFilterParams.TauDamp = 0.03;
surroundFilterParams.Period = 0.5;
surroundFilterParams.Phase = 180; % % %
surroundFilterParams.ScaleFactor = 10;

centerNLParams.alphaScale = 0.5;
centerNLParams.betaSens = 4;
centerNLParams.gammaXoffset = -0.4;
centerNLParams.epsilonYoffset = -0.1;

surroundNLParams.alphaScale = 0.5;
surroundNLParams.betaSens = 2;
surroundNLParams.gammaXoffset = -0.4;
surroundNLParams.epsilonYoffset = -0.4;

stimTime = 0.001:0.001:1;
centerStimulus = zeros(size(stimTime));
surroundStimulus = zeros(size(stimTime));

centerStimulus(50:150) = 1;
surroundStimulus(150:250) = 1;

centerStimulus(300:400) = 1;
surroundStimulus(300:400) = 1;

centerStimulus(450:550) = -1;
surroundStimulus(450:550) = -1;

centerStimulus(600:700) = 1;
surroundStimulus(600:700) = -1;

centerStimulus(750:850) = -1;
surroundStimulus(750:850) = 1;


[response, modelComponents] = ParamCSModel_independent(centerStimulus,surroundStimulus,...
    centerFilterParams,surroundFilterParams,centerNLParams,surroundNLParams);
response = response + 0.1.*randn(size(response));

figure(1); clf; subplot(3,1,1); plot(centerStimulus,'b');
subplot(312); plot(surroundStimulus,'r')
subplot(313); plot(response);


%%
% % coef = [10, 0.1, 0.03, 0.5, -10, 10,...
% %     10, 0.1, 0.03, 0.5, 180, 10,...
% %     0.5, 4, -0.4, -0.1,...
% %     0.5, 4, -0.4, -0.1];

coef = [1, 0.4, 0.01, 0.5, -20, 10,...
    20, 0.05, 0.01, 0.5, 100, 10,...
    0.1, 3, -0.2, -0.1,...
    0.4, 1, -0.4, -0.8];


for iter = 1:4
    fitcoef = IndepModelWrapper(coef, [centerStimulus; surroundStimulus], response);
    coef = fitcoef;
end
