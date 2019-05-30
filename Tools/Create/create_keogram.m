function [keo,parallels,meridian] = create_keogram(image,latitude,longitude,varargin)
%create_keogram.m Creates the keogram, at the middle-meridian or the
%meridian you specify, and also returns the corresponding parallels.

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter(p,'meridian',nan);
addParameter(p,'latPixelNum',nan);
addParameter(p,'projectionAlt',nan);
addParameter(p,'time',nan);
addParameter(p,'sensorLoc',nan); % [lat, lon, alt]

addRequired(p,'image');
addRequired(p,'latitude',@(x)ismatrix(x));
addRequired(p,'longitude',@(x)ismatrix(x));


parse(p,image,latitude,longitude,varargin{:});

meridian = p.Results.meridian;



nTime = size(p.Results.image,1);
image = reshape(p.Results.image,nTime,[])';
lat   = reshape(p.Results.latitude,1,[])';
lon   = reshape(p.Results.longitude,1,[])';

latPixelNum = p.Results.latPixelNum; 

if isnan(latPixelNum) % Choosing number of latitude pixels
    latPixelNum = ceil(sqrt(length(lat)));
end

parallels = linspace(min(lat),max(lat),latPixelNum);

if isnan(meridian) % Chosing the keogram meridian to be the magnetic meridian
    if isnan(p.Results.sensorLoc)
        sensorLoc(1) = 65.126;
        sensorLoc(2) = -147.4789;
        sensorLoc(3) = 0.689;
    else
        sensorLoc = p.Results.sensorLoc;
    end
    
    if isnan(p.Results.time)
        time = datenum('26 Mar 2008 11:00');
    else
        time = p.Results.time;
    end
    
    if isnan(p.Results.projectionAlt)
        projectionAlt = 110;
    else
        projectionAlt = p.Results.projectionAlt;
    end
%     meridian = (min(longitude(:))+max(longitude(:)))/2;
    meridian = get_magnetic_meridian(sensorLoc,time,parallels,projectionAlt);
elseif numel(meridian)==1
    meridian = meridian.*ones(size(parallels));
end

multiWaitbar('Keogramming...',0);  
id = 1./nTime;
    for i=1:1:nTime
        F = scatteredInterpolant(lat(~isnan(lat)),lon(~isnan(lat)),...
            image(~isnan(lat),i));
        keo(:,i) = F(parallels,meridian);
        multiWaitbar('Keogramming...','Increment',id);
    end
multiWaitbar('Keogramming...',1);  
end


