function p = create_panels(figureHandle,varargin)
%create_panels.m Creates a set of panels on a Letter Size page
%
% Input
%--------
% figureHandle - input figure handle
% totalPanelNo - Total number of panels
% panelHeight    - Height of the panel in [mm], Default: 30
% demargin     - Margin between each panel in [mm], Default: 4
% marginleft   - Left margin in [mm], Default: 35
% marginright  - Right margin 
% margintop    - Top margin 

po = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);


addParameter(po,'totalPanelNo',1,validScalarPosNum); % in mmm
addParameter(po,'panelHeight',30,validScalarPosNum); % in mm
addParameter(po,'panelBreadth',[],validScalarPosNum); % in mm
addParameter(po,'demargin',4,validScalarPosNum); % in mm
addParameter(po,'marginleft',35,validScalarPosNum);
addParameter(po,'marginright',25,validScalarPosNum);
addParameter(po,'margintop',10,validScalarPosNum);
addRequired(po,'figureHandle',@(x)isfigure(x));

parse(po,figureHandle,varargin{:});

% clf
p=panel();
p.pack(1);

panelHeight = po.Results.panelHeight; %in mm
demargin = po.Results.demargin;
totalPanelNo = po.Results.totalPanelNo;

panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelHeight}}; 
end

p(1).pack(panelDefinition);

p.marginleft=po.Results.marginleft;
p.marginright=po.Results.marginright;
p.margintop = po.Results.margintop;
p(1).de.margin=demargin;
% p.fontsize=12;
if ~isempty(po.Results.panelBreadth)
    panelBreadth = po.Results.panelBreadth+po.Results.marginleft+po.Results.marginright;
else
    panelBreadth = [];
end
resize_figure(figureHandle, (panelHeight+demargin)*totalPanelNo+4*(totalPanelNo+1),...
    panelBreadth);
p.select('all');
  


end

function OK = isfigure(h)
 if strcmp(get(h,'type'), 'figure')
     OK = true;
 else
     OK = false;
 end
     
end