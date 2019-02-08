function [meridian] = get_magnetic_meridian(sensorLoc,time,parallel,projectionAltitude)
% Input
% sensorLoc - [lat(deg), lon(deg), alt(km)
% time      - matlab time
% parallel  - an array of latitude values in degrees, for which you need
%             the longitude values of the magnetic meridian
% meridian  - the meridian values of the magnetic meridian for the input
%             parallels
% projectionAltitude - km;

[Bx, By] = igrf(time,sensorLoc(1),sensorLoc(2),projectionAltitude);
% calculating declination: the azimuth from true north (East +ve)
declination = rad2deg(atan2(By,Bx));

elevation = 0:0.1:90;

el = [elevation(1:end-1),elevation(end:-1:1)];

az= [repmat(declination+180,1,length(elevation)-1),repmat(declination,1,length(elevation))];

RE = 6.371*10^3;
projectionAltitude1 = projectionAltitude - sensorLoc(3);
slantRange = -RE.*sind(el) + modSign(el).*sqrt((RE^2).*(sind(el)).^2 + projectionAltitude1.*(projectionAltitude1+2.*RE)); 
slantRange(slantRange<=projectionAltitude1*0.1) = nan;

[lat,lon,alt] = aer2geodetic(az,el,slantRange,sensorLoc(1),sensorLoc(2),sensorLoc(3),wgs84Ellipsoid('km'));

meridian = interp1(lat, lon, parallel, 'linear');

end

function y = modSign(x)
    for i=1:1:length(x)
        if x(i)==0
            x(i) = +1;
        end
    end
    y = sign(x);
end