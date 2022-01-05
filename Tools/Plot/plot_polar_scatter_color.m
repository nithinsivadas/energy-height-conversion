function [pax] = plot_polar_scatter(mlat,mlt,varargin)
%plot_polar_scatter Plots a polar scatter plot with input mlat and mlt 
%   Detailed explanation goes here

%% Parse Inputs
p = inputParser;
addParameter(p,'RLim',[0,90]);

addRequired(p,'mlat');
addRequired(p,'mlt');

parse(p,mlat,mlt,varargin{:});

%% Plot
polarscatter(deg2rad(mlt.*360/24),mlat,1,'filled');

%% Set axis properties
pax = gca;
pax.ThetaTickLabel=num2cell(pax.ThetaTick.*24/360);
pax.RDir='reverse';
pax.ThetaDir = 'counterclockwise';
pax.ThetaZeroLocation = 'bottom';
pax.RLim = p.Results.RLim;
end

