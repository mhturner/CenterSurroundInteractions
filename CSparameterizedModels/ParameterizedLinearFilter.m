function response = ParameterizedLinearFilter(timeVector, varargin)
    ip = inputParser;
    ip.addRequired('timeVector',@isnumeric); %sec
    
    addParameter(ip,'NumFilters',10,@isnumeric);
    addParameter(ip,'TauRise',0.1,@isnumeric); %sec
    addParameter(ip,'TauDamp',0.03,@isnumeric); %sec
    addParameter(ip,'Period',0.5,@isnumeric);
    addParameter(ip,'Phase',-10,@isnumeric); %degrees
    addParameter(ip,'ScaleFactor',10,@isnumeric);
    
    ip.parse(timeVector,varargin{:});
    timeVector = ip.Results.timeVector;
    
    NumFilters = ip.Results.NumFilters;
    TauRise = ip.Results.TauRise;
    TauDamp = ip.Results.TauDamp;
    Period = ip.Results.Period;
    Phase = ip.Results.Phase;
    ScaleFactor = ip.Results.ScaleFactor;

    response = ScaleFactor .* (((timeVector./TauRise).^NumFilters)./(1+((timeVector./TauRise).^NumFilters))) .* ...
        exp(-((timeVector./TauDamp))) .* cos(((2.*pi.*timeVector)./Period)+(2*pi*Phase/360));
end