function [keo,parallels,meridian] = create_keogram(image,latitude,longitude,varargin)
%create_keogram.m Creates the keogram, at the middle-meridian or the
%meridian you specify, and also returns the corresponding parallels.

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter(p,'meridian',nan);
addParameter(p,'latPixelNum',nan);

addRequired(p,'image');
addRequired(p,'latitude',@(x)ismatrix(x));
addRequired(p,'longitude',@(x)ismatrix(x));

parse(p,image,latitude,longitude,varargin{:});

meridian = p.Results.meridian;

if isnan(meridian) % Chosing the keogram meridian to be the center longitude
    meridian = (min(longitude(:))+max(longitude(:)))/2;
end

nTime = size(p.Results.image,1);
image = reshape(p.Results.image,nTime,[])';
lat   = reshape(p.Results.latitude,1,[])';
lon   = reshape(p.Results.longitude,1,[])';

latPixelNum = p.Results.latPixelNum; 

if isnan(latPixelNum) % Choosing the keogram meridian to be the center longitude
    latPixelNum = ceil(sqrt(length(lat)));
end

parallels = linspace(min(lat),max(lat),latPixelNum);
multiWaitbar('Keogramming...',0);  
id = 1./nTime;
    for i=1:1:nTime
        F = scatteredInterpolant(lat(~isnan(lat)),lon(~isnan(lat)),...
            image(~isnan(lat),i));
        keo(:,i) = F(parallels,meridian*ones(size(parallels)));
        multiWaitbar('Keogramming...','Increment',id);
    end
multiWaitbar('Keogramming...',1);  
end

