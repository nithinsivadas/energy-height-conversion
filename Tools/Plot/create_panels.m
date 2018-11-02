function p = create_panels(figureHandle,varargin)
%create_panels.m Creates a set of panels on a Letter Size page
%
% Input
%--------
% figureHandle - input figure handle
% totalPanelNo - Total number of panels
% panelSize    - Height of the panel in [mm], Default: 30
% demargin     - Margin between each panel in [mm], Default: 4
% marginleft   - Left margin in [mm], Default: 35
% marginright  - Right margin 
% margintop    - Top margin 

po = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);


addParameter(po,'totalPanelNo',1,validScalarPosNum); % in mmm
addParameter(po,'panelSize',30,validScalarPosNum); % in mmm
addParameter(po,'demargin',4,validScalarPosNum); % in mmm
addParameter(po,'marginleft',35,validScalarPosNum);
addParameter(po,'marginright',25,validScalarPosNum);
addParameter(po,'margintop',10,validScalarPosNum);
addRequired(po,'figureHandle',@(x)isfigure(x));

parse(po,figureHandle,varargin{:});

% clf
p=panel();
p.pack(1);

panelSize = po.Results.panelSize; %in mm
demargin = po.Results.demargin;
totalPanelNo = po.Results.totalPanelNo;

panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelSize}}; 
end

p(1).pack(panelDefinition);

p.marginleft=po.Results.marginleft;
p.marginright=po.Results.marginright;
p.margintop = po.Results.margintop;
p(1).de.margin=demargin;
% p.fontsize=12;
resize_figure(figureHandle, panelSize*totalPanelNo+4*(totalPanelNo+1)+20);
p.select('all');
  


end

function OK = isfigure(h)
 if strcmp(get(h,'type'), 'figure')
     OK = true;
 else
     OK = false;
 end
     
end