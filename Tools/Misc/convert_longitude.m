function [longitude] = convert_longitude(longitude,modeStr)
%convert_longitude Converts longitude from 0-360 to -180(W)to+180(E), and vice
%versa. 
%   Detailed explanation goes here
if nargin<2
    modeStr = '360to180';
end
flag = 0;
if strcmp(modeStr,'360to180')
    if sum(longitude<0)<1
        longitude = rem((longitude+180),360)-180;
    else
        flag = 1;
    end
elseif strcmp(modeStr,'180to360')
    if sum(abs(longitude)>180)<1
        longitude =  mod(longitude,360);
    else
        flag = 1;
    end
else
    flag = 1;
end

if flag == 1
    error('The modestr is not appropriate for the data');
end

end

