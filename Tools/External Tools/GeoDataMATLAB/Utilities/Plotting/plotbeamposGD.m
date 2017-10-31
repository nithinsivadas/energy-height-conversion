function plotbeamposGD(GD,varargin)
% plotbeamposGD.m
% by John Swoboda
% This function will plot the beam grid in a polar grid with the az angle
% as the azimithal commponent of the grid and the el angle as the range in
% decending order.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% GD - An instance of the geodata class in Spherical coordinates.
%
% Properties - These are Name Value pairs that are optional inputs like in 
% MATLAB's plotting function. The name will be listed first and the
% internal variable name will be stated in the description.
% 
% Fig- The value will be a MATLAB figure handle. The internal variable name
% is figname. The default value is nan
% axh- The value will be a MATLAB axis handle. The internal variable name
% is axh. The default value is nan
% title- The value will be a MATLAB string for the plot title. The internal
% variable name is titlestr. The default value is an empty string.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
% Determine the properties
paramstr = varargin(1:2:end);
paramvals = varargin(2:2:end);
poss_labels={'Fig','axh','title'};
varnames = {'figname','axh','titlestr'};
vals = {nan,nan,'Beampositions'};
checkinputs(paramstr,paramvals,poss_labels,vals,varnames)


if isnumeric(figname)
    figure();
    axh = newplot(figname);
else
    figure(figname);
end
if isnumeric(axh);
    axh=gca;
end

%% Determine angles

assert(strcmpi(GD.coordnames,'spherical'),'The GeoData instance must be in spherical coordinates');
angles = GD.dataloc(:,2:end);
uang = unique(angles,'rows');

az = uang(:,1);
el = uang(:,2);
az = az*pi/180;


%% Plotting

theta = linspace(0,2*pi,100);
r10 = 70*ones(1,100);
[xx1,yy1] = pol2cart(theta,r10);
fill(xx1,yy1,'w')
hold on

axis off 
axis equal
hold on

%plotting the radar beams
[xx,yy]=pol2cart(-az+pi/2,90-el);
plot(xx,yy,'o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',10)


%Plotting the cardinal directions


lin1x=linspace(-70,70,700);
lin1y=linspace(0,0,700);
plot(lin1x,lin1y,':k')
lin1x=linspace(0,0,700);
lin1y=linspace(-70,70,700);
plot(lin1x,lin1y,':k')
lin1x=linspace(-70*sind(45),70*sind(45),700);
lin1y=linspace(-70*cosd(45),70*cosd(45),700);
plot(lin1x,lin1y,':k')
lin1x=linspace(-70*sind(45),70*sind(45),700);
lin1y=linspace(70*cosd(45),-70*cosd(45),700);
plot(lin1x,lin1y,':k')

text(73,0,'E')
text(73*cosd(45),73*sind(45),'NE')
text(-1,73,'N')
text(-76,0,'W')
text(-81*cosd(45),76*sind(45),'NW')
text(-1,-73,'S')
text(-81*cosd(45),-76*sind(45),'SW')
text(73*cosd(45),-73*sind(45),'SE')

%plotting elevation rings
elev=[10,20,30,40,50,60];
theta = linspace(0,2*pi,100);
r20 = elev'*ones(1,100);
for i=1:length(elev)
    [x,y] = pol2cart(theta,r20(i,:));
    plot(x,y,':k')
    text(elev(i)*sind(22.5),-elev(i)*cosd(22.5)-2,sprintf('%d',90-elev(i)))
end
title({titlestr;''},'FontSize',16);
axis off
% 
% for i=1:1:length(xx)
%     text(xx(i)+2,yy(i)+2,sprintf('%d',i))     
% end       
        