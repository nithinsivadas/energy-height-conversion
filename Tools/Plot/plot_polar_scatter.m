function [pax] = plot_polar_scatter(mlat,mlt,varargin)
%plot_polar_scatter Plots a polar scatter plot with input mlat and mlt 
%   Detailed explanation goes here

%% Parse Inputs
p = inputParser;
addParameter(p,'RLim',[0,90]);
addParameter(p,'zColor',nan);
addParameter(p,'zLim',nan);
addParameter(p,'markerSize',5);
addParameter(p,'colorMin',nan);
addParameter(p,'colorMax',nan);

addRequired(p,'mlat');
addRequired(p,'mlt');

parse(p,mlat,mlt,varargin{:});

if isnan(p.Results.zColor)
    %% Plot
    polarscatter(deg2rad(mlt.*360/24),mlat,1,'filled');

    %% Set axis properties
    pax = gca;
    pax.ThetaTickLabel=num2cell(pax.ThetaTick.*24/360);
    pax.RDir='reverse';
    pax.ThetaDir = 'counterclockwise';
    pax.ThetaZeroLocation = 'bottom';
    pax.RLim = p.Results.RLim;
else
    %% Plot
    zData = p.Results.zColor;
    if isnan(p.Results.zLim)
        zLim = [min(zData), max(zData)];
    else
        zLim = p.Results.zLim;
    end
    zN = length(zData);
    zLinear = linspace(zLim(1),zLim(2),zN)';
    cmap = get_colormap3('b','g','r',zN);
    cData = interp1(zLinear,cmap,zData);
    colormap(cmap);
    polarscatter(deg2rad(mlt.*360/24),mlat,p.Results.markerSize,cData,'filled','Marker','o');
    %% Set axis properties
    pax = gca;
    pax.ThetaTickLabel=num2cell(pax.ThetaTick.*24/360);
    pax.RDir='reverse';
    pax.ThetaDir = 'counterclockwise';
    pax.ThetaZeroLocation = 'bottom';
    pax.RLim = p.Results.RLim;
    colorbar;
    caxis([zLinear(1) zLinear(end)]);
end
end

