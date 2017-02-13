function plot_DASC_aer( dataNew, az, el, imageSize )
%% plot_DASC_aer.m Plots the All Sky Image in az, el, range coordinates
%--------------------------------------------------------------------------
% Input
%------
% dataNew - 2-D optical data from DASC in geodetic coordinates [nCoordinates]
% az     - Azimuth coordinates [nCoordinates]
% el     - Elevation coordinates [nCoordinates]
% imageSize - the total pixel size of the output image e.g. 1024 
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
% include_new_colormaps;
% colormap(viridis);

% Rotating az data 90 deg counter-clockwise [Very Important]
x = (90-el).*sind(rotate_array(az(:),-90));
y = (90-el).*cosd(rotate_array(az(:),-90)); 


F1 = scatteredInterpolant(x,y,dataNew,'linear','none');

xq = linspace(-90,90,imageSize);
yq = linspace(-90,90,imageSize);
[X,Y] =ndgrid(xq,yq);

Vq1 = F1(X, Y);

h2=pcolor(X,Y,Vq1);
set(h2,'EdgeColor','none');
colorbar;
axis equal
% hold on; annotation('textarrow',[0.15 0.20],[0.8 0.8],'string','W');
% hold on; annotation('textarrow',[0.2 0.2],[0.85 0.90],'string','N');
hold on; text(-77,75,'N');
hold on; line([-75 -75],[70 60],'color','r');
hold on; text(-77,55,'S');

hold on; text(-69,65,'W');
hold on; line([-80 -70],[65 65],'color','r');
hold on; text(-85,65,'E');

end

