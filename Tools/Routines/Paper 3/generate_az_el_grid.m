
imageSize = 1024;

x = linspace(-90,90,imageSize);
y = linspace(-90,90,imageSize);

[X, Y] = meshgrid(x,y);

az = wrapTo360(atan2d(X,Y));
el = 90 - X./(sind(az));

indx = el<=0;
az(indx) = nan;
el(indx) = nan;