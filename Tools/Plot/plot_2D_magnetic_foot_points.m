function [ax2, h2] = plot_2D_magnetic_foot_points...
    (magEqPointGEO, ionosphereCoord, BfieldModelStr, thisTimeBfieldModel, setMapOn, latLim, lonLim,...
    contourArray,setFieldLabelOn, setTimeLabelOn, imageSize)
%plot_2D_magnetic_foot_points Given a matrix of magnetic conjugate points
%at the equatorial plane, and their corresponding ionospheric foot points,
%the function draws a contour plot of their equatorial radial distances.
%   Detailed explanation goes here

if nargin<11 || isempty(imageSize)
    imageSize = 512; end

if nargin < 10 || isempty(setTimeLableOn)
        setTimeLabelOn = true; end

if nargin < 9 || isempty(setFieldLabelOn)
    setFieldLabelOn = true; end

if nargin < 8 || isempty(contourArray)
    contourArray = 1:1:30; end

if nargin < 5 ||isempty(setMapOn)
    setMapOn = true; end

if nargin < 4 ||isempty(thisTimeBfieldModel)
thisTimeBfieldModel = datenum('9 September 9999'); end

if nargin < 3 ||isempty(BfieldModelStr)
BfieldModelStr = 'Unknown'; end

RE=(magEqPointGEO(1,:).^2+magEqPointGEO(2,:).^2+magEqPointGEO(3,:).^2).^0.5;
GDZ = ionosphereCoord;

F = scatteredInterpolant(GDZ(:,2),GDZ(:,3),RE');

if nargin < 7 ||isempty(lonLim)
    lonLim = [min(GDZ(:,3)) max(GDZ(:,3))]; end

if nargin < 6 ||isempty(latLim)
    latLim = [min(GDZ(:,2)) max(GDZ(:,2))]; end

imageSize = 512;
latq = linspace(latLim(1),latLim(2),imageSize);
lonq = linspace(lonLim(1),lonLim(2),imageSize);

Vq = F({latq, lonq});
Vq(Vq<=0)=nan;
    
    latWidth = latLim(2) - latLim(1);
    lonWidth = lonLim(2) - lonLim(1);
    latCenter = mean(latLim);
    lonCenter = mean(lonLim);
    
if setMapOn==true
    ax2=axesm('lambertstd','MapLatLimit',[(latLim(1)) (latLim(2))],...
        'MapLonLimit',[(lonLim(1)) (lonLim(2))],...
        'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
        'PLineLocation',2,'MLineLocation',5);
%     ax2=axesm('lambertstd','FLatLimit',[(latLim(1))-latWidth/2 (latLim(2))+latWidth/2],...
%         'FLonLimit',[(lonLim(1))-lonWidth/2 (lonLim(2))+lonWidth/2],...
%         'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on',...
%         'PLineLocation',2,'MLineLocation',5);    
    axis off
    load coastlines
    plotm(coastlat,coastlon)
    hold on;
end
ax = gca;
xlim = ax.XLim;
ylim = ax.YLim;
[c2,h2]=contourm(latq,lonq,Vq,contourArray,'LineColor',[0 1 1]); 

if setFieldLabelOn == true
    clabelm(c2,h2,'LabelSpacing',723);
end

ax.XLim = xlim;
ax.YLim = ylim;

if setTimeLabelOn==true
    hold on;
    textm(latLim(2)-0.75*latWidth, lonLim(2)-0.3*lonWidth, BfieldModelStr,'color','c');
    textm(latLim(2)-0.75*latWidth+latWidth/20, lonLim(2)-0.3*lonWidth,...
        [datestr(thisTimeBfieldModel,'HH:MM:SS'),' UT'],'color', 'c');
    hold off;
end

end

