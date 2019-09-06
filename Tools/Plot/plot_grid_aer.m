function plot_grid_aer(az, el, colorStr)
%% plot_grid_aer.m Plots grid lines in polar plot at specified azimuth and elevation
%--------------------------------------------------------------------------
% Input
%------
% az        - Array of azimuth values to plot grid lines
% el        - Array of elevation values to plot grid lines
% colorStr  - Color of the grid lines e.g. 'r' 
%--------------------------------------------------------------------------
% Modified: 24th Jan 2017 
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%--------------------------------------------------------------------------
az_all = 0:1:360;
el_all = 0:1:90;

if ~isnan(el)
    for iel = 1:1:length(el)
        x = (90-el(iel)).*sind(360-az_all); %Counter clockwise rotation
        y = (90-el(iel)).*cosd(az_all);
        hold on;
        plot(x,y,colorStr);
        text(max(x),10,[num2str(el(iel)),'^0'],'color',colorStr);
    end
end

if ~isnan(az)
%      az = rotate_array(az,90); % Since 0-> west, making it point to north
    for iaz = 1:1:length(az)
        x = (90-el_all).*sind(360-az(iaz));
        x1 = (90-el_all).*sind(rotate_array(360-az(iaz),180));
        y = (90-el_all).*cosd(az(iaz));
        y1 = (90-el_all).*cosd(rotate_array(az(iaz),180));
        hold on;
        plot(x,y,colorStr);
        plot(x1,y1,colorStr);
        [mx,i]=max(abs(x));
        [my,j]=max(abs(y));
        text((mx+5)*sign_user(x(i)),(my+6)*sign_user(y(j)),[num2str(rotate_array(az(iaz),0)),'^0'],'color',colorStr);
        [mx1,i]=max(abs(x1));
        [my1,j]=max(abs(y1));
        text((mx1-2)*sign_user(x1(i)),(my1+6)*sign_user(y1(j)),[num2str(rotate_array(rotate_array(az(iaz),0),180)),'^0'],'color',colorStr);
    end
end    
end

function y= sign_user(x)
    if x~=0
        y = sign(x);
    else
        y = 1;
    end
end