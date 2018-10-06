function [ax2, h2, c2] = plot_2D_magnetic_foot_points...
    (plotVariable, ionosphereCoord, varargin)
%plot_2D_magnetic_foot_points Given a matrix of magnetic conjugate points
%at the equatorial plane, and their corresponding ionospheric foot points,
%the function draws a contour plot of their equatorial radial distances.
%   Detailed explanation goes here
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

addParameter(p,'imageSize', 512, validScalarPosNum);
addParameter(p,'setTimeLabelOn', true);
addParameter(p,'setFieldLabelOn', true);
addParameter(p,'contourLineArray', 1:1:30);
addParameter(p,'contourLabelArray', 1:4:30);
addParameter(p,'latLim', 0);
addParameter(p,'lonLim', 0);
addParameter(p,'setMapOn', true);
addParameter(p,'thisTimeBfieldModel', datenum('9 September 9999'));
addParameter(p,'BfieldModelStr','Unknown');

addRequired(p,'plotVariable'); % Can be L-shell, L*, or anything else
addRequired(p,'ionosphereCoord');

parse(p,plotVariable, ionosphereCoord, varargin{:});

GDZ = ionosphereCoord;

F = scatteredInterpolant(GDZ(:,2),GDZ(:,3),plotVariable);

if p.Results.lonLim == 0
    lonLim = [min(GDZ(:,3)) max(GDZ(:,3))];
else
    lonLim = p.Results.lonLim;
end

if p.Results.latLim == 0
    latLim = [min(GDZ(:,2)) max(GDZ(:,2))]; 
else
    latLim = p.Results.latLim;
end

latq = linspace(latLim(1),latLim(2),p.Results.imageSize);
lonq = linspace(lonLim(1),lonLim(2),p.Results.imageSize);

Vq = F({latq, lonq});
Vq(Vq<=0)=nan;
    
    latWidth = latLim(2) - latLim(1);
    lonWidth = lonLim(2) - lonLim(1);
    latCenter = mean(latLim);
    lonCenter = mean(lonLim);
    
if p.Results.setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1)) (latLim(2))],...
        'MapLonLimit',[(lonLim(1)) (lonLim(2))],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',2,'MLineLocation',5);
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
end
ax = gca;
xlim = ax.XLim;
ylim = ax.YLim;
[c2,h2]=contourm(latq,lonq,Vq,p.Results.contourLineArray,'LineColor',[0 1 1]); 

if p.Results.setFieldLabelOn == true
    htext = clabelm(c2,h2,p.Results.contourLabelArray,...
        'LabelSpacing',723);
    set(htext,'BackgroundColor','cyan','margin',2);
end

ax.XLim = xlim;
ax.YLim = ylim;

if p.Results.setTimeLabelOn==true
    hold on;
    textm(latLim(2)-0.75*latWidth, lonLim(2)-0.3*lonWidth, [p.Results.BfieldModelStr,' L-shell'],'color','c');
    textm(latLim(2)-0.75*latWidth+latWidth/20, lonLim(2)-0.3*lonWidth,...
        [datestr(p.Results.thisTimeBfieldModel,'HH:MM:SS'),' UT'],'color', 'c');
    hold off;
end

end

