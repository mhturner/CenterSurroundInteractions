function response = ParameterizedNonlinearity(Input,varargin)
    ip = inputParser;
    ip.addRequired('Input',@isnumeric);
    
    addParameter(ip,'alphaScale',0.5,@isnumeric);
    addParameter(ip,'betaSens',4,@isnumeric);
    addParameter(ip,'gammaXoffset',-0.4,@isnumeric);
    addParameter(ip,'epsilonYoffset',-0.1,@isnumeric);
    
    ip.parse(Input,varargin{:});
    Input = ip.Results.Input;
    
    alphaScale = ip.Results.alphaScale;
    betaSens = ip.Results.betaSens;
    gammaXoffset = ip.Results.gammaXoffset;
    epsilonYoffset = ip.Results.epsilonYoffset;

    response = alphaScale*normcdf(betaSens.*Input + gammaXoffset,0,1)+epsilonYoffset;
end