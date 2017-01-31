function fitcoef = IndepModelWrapper(coef, stimulusMatrix, response)
    fitOptions = optimset('MaxIter',1500,'MaxFunEvals',5000);
    fitcoef = lsqnonlin(@IndepModelErr, coef, [], [], fitOptions, stimulusMatrix, response);
end

function err = IndepModelErr(coef,stimulusMatrix,response)
    centerFilterParams.NumFilters = coef(1);
    centerFilterParams.TauRise = coef(2);
    centerFilterParams.TauDamp = coef(3);
    centerFilterParams.Period = coef(4);
    centerFilterParams.Phase = coef(5);
    centerFilterParams.ScaleFactor = coef(6);

    surroundFilterParams.NumFilters = coef(7);
    surroundFilterParams.TauRise = coef(8);
    surroundFilterParams.TauDamp = coef(9);
    surroundFilterParams.Period = coef(10);
    surroundFilterParams.Phase = coef(11);
    surroundFilterParams.ScaleFactor = coef(12);

    centerNLParams.alphaScale = coef(13);
    centerNLParams.betaSens = coef(14);
    centerNLParams.gammaXoffset = coef(15);
    centerNLParams.epsilonYoffset = coef(16);

    surroundNLParams.alphaScale = coef(17);
    surroundNLParams.betaSens = coef(18);
    surroundNLParams.gammaXoffset = coef(19);
    surroundNLParams.epsilonYoffset = coef(20);
    
    centerStimulus = stimulusMatrix(1,:);
    surroundStimulus = stimulusMatrix(2,:);

    [fitResponse, ~ ] = ParamCSModel_independent(centerStimulus,surroundStimulus,centerFilterParams,surroundFilterParams,centerNLParams,surroundNLParams);

    err = sum((response - fitResponse).^2);
end