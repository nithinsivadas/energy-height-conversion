function [obslat, obslong, obsalt] = getobs

% interactive request of observer coordinates

% output

%  obslat  = geographic latitude (radians)
%  obslong = geographic longitude (radians)
%  obsalt  = geodetic altitude (kilometers)

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180;

for itry = 1:1:5
    fprintf('\nplease input the geographic latitude of the observer');
    fprintf('\n(-90 <= degrees <= +90, 0 <= minutes <= 60, 0 <= seconds <= 60)');
    fprintf('\n(north latitude is positive, south latitude is negative)\n');

    latstr = input('? ', 's');

    tl = size(latstr);

    ci = findstr(latstr, ',');

    % extract degrees, minutes and seconds

    latdeg = str2double(latstr(1:ci(1)-1));

    latmin = str2double(latstr(ci(1)+1:ci(2)-1));

    latsec = str2double(latstr(ci(2)+1:tl(2)));

    % check for valid inputs

    if (abs(latdeg) >= 0 && abs(latdeg) <= 90 && ...
            latmin >= 0 && latmin <= 60 && ...
            latsec >= 0 && latsec <= 60)
        break;
    end
end

for itry = 1:1:5
    fprintf('\nplease input the geographic longitude of the observer');
    fprintf('\n(0 <= degrees <= 360, 0 <= minutes <= 60, 0 <= seconds <= 60)');
    fprintf('\n(east longitude is positive, west longitude is negative)\n');

    longstr = input('? ', 's');

    tl = size(longstr);

    ci = findstr(longstr, ',');

    % extract degrees, minutes and seconds

    longdeg = str2double(longstr(1:ci(1)-1));

    longmin = str2double(longstr(ci(1)+1:ci(2)-1));

    longsec = str2double(longstr(ci(2)+1:tl(2)));

    % check for valid inputs

    if (abs(longdeg) >= 0 && abs(longdeg) <= 360 && ...
            longmin >= 0 && longmin <= 60 && ...
            longsec >= 0 && longsec <= 60)
        break;
    end
end

fprintf('\nplease input the altitude of the observer (meters)');
fprintf('\n(positive above sea level, negative below sea level)\n');

obsalt = input('? ');

% convert latitude to radians

if (latstr(1:1) == '-' && latdeg == 0)
    obslat = -dtr * (latmin / 60 + latsec / 3600);
elseif (latdeg == 0)
    obslat = dtr * (latmin / 60 + latsec / 3600);
else
    obslat = dtr * sign(latdeg) * (abs(latdeg) + latmin / 60 ...
        + latsec / 3600);
end

% convert longitude to radians

if (longstr(1:1) == '-' && longdeg == 0)
    obslong = -dtr * (longmin / 60 + longsec / 3600);
elseif (longdeg == 0)
    obslong = dtr * (longmin / 60 + longsec / 3600);
else
    obslong = dtr * sign(longdeg) * (abs(longdeg) + longmin / 60 ...
        + longsec / 3600);
end

% convert altitude to kilometers

obsalt = 0.001 * obsalt;
