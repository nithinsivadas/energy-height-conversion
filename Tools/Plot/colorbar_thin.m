function [c] = colorbar_thin(varargin)
%colorbar_thin.m Creates thin colorbar
%   Detailed explanation goes here
p = inputParser;
addParameter(p,'Location','eastoutside');
addParameter(p,'YLabel','');
addParameter(p,'Width',0.2);

parse(p,varargin{:});

axPos = get(gca, 'position');
c = colorbar(p.Results.Location);
set(gca,'position',axPos);

cPos=get(c,'Position');
cPos(3)=p.Results.Width*cPos(3);
set(c, 'Position',cPos);

ylabel(c,p.Results.YLabel,'Interpreter','Latex');

end

