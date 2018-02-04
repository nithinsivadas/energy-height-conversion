function hslice = sliceGD(GD,varargin)
% sliceGD.m
% by John Swoboda
% This function will access the three-d slice function for the GeoData
% objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% GD - An instance of the geodata class in Cartisian coordinates.
% sx,sy,sz - Vectors where volumetric slices will be drawn.
% Xi,Yi,Zi - Matrices that hold thr surfaces that that data will be plotted
% over.
% Properties - These are Name Value pairs that are optional inputs like in 
% MATLAB's plotting function. The name will be listed first and the
% internal variable name will be stated in the description.

% key - A string of the variable name in the data set that will be plotted.
% The internal variable name is key and the default value will be the first 
% string in the data struct of GD.
% Fig- The value will be a MATLAB figure handle. The internal variable name
% is figname. The default value is nan and new figure will be created.
% axh- The value will be a MATLAB axis handle. The internal variable name
% is axh. The default value is nan and a new axis will be created.
% title- The value will be a MATLAB string for the plot title. The internal
% variable name is titlestr. The default value is an empty string.
% time - The value will be an integer that cooresponds to the time value in
% GD. The internal variable name is timenum and the default value is 1.
% bounds - A vector of two elements, the first is the lower bound and the
% second is the higher bound for the Caxis. The internal variable name is 
% vbound and its default value will be the max and min of the data to be
% plotted.
% twodplot - This is a bool that will determine if the data is two
% dimensional and needs to be plotted over a surface. The internal variable
% name is twodplot and its default value is false.
% colormap - This the colormap for the data. The internal variable name is
% cmap and its default value is MATLAB's default colormap.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% hslice - The handle for the plotted object.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
assert(strcmpi(GD.coordnames,'cartesian'),'GeoData object needs to be in cartisian coordinates.')
% Determine if given a surface or specific slice values.
emptyarr = false(1,3);
for k=1:3
    emptyarr(k)= isempty(varargin{k});
end
if (ndims(varargin{1})==1)||any(emptyarr)
    [sx,sy,sz] = varargin{1:3};
    surftype = 2;
elseif ismatrix(varargin{1})
    [Xi,Yi,Zi] = varargin{1:3};
    surftype = 1;

end  

% determine which is the default colormap
v2014dt = datetime('September 15, 2014');
[~,d] = version();

if datetime(d)>=v2014dt
    defmap = parula(64);
else
    defmap = jet(64);
end
% Determine the properties
paramstr = varargin(4:2:end);
paramvals = varargin(5:2:end);
poss_labels={'key','Fig','axh','title','time','bounds','twodplot','colormap'};
varnames = {'key','figname','axh','titlestr','timenum','vbound','twodplot','cmap'};
vals = {1,nan,nan,'Generic',1,[nan,nan],false,defmap};
checkinputs(paramstr,paramvals,poss_labels,vals,varnames)

% apply default parameters 
if isnumeric(key)
    dnames = fieldnames(GD.data);
    key=dnames{key};
end

if isnumeric(figname)
    figname = figure();
    axh = newplot(figname);
else
    figure(figname);
end
if isnumeric(axh);
    axh=gca;
end

    
% Augment the title string to remove the wildcard characters
titlestr = insertinfo(titlestr,'key',key,'time',GD.times(timenum,1),'timend',GD.times(timenum,2));

v = GD.data.(key)(:,timenum);
if all(isnan(vbound))
    vbound = [min(v),max(v)];
end
% for a twod object apply a colormap to the surface because the slice
% function requires a volumetric object
if twodplot
    if surftype==2
        % reshape the data
        [X,Y,Z,V] = reshapegen(GD.dataloc,v);
        origloc = GD.dataloc(1,:);
        ONlocs = size(GD.dataloc,1);
        loclog = origloc(ones(ONlocs,1),:)==GD.dataloc;
        %XXX Need to come up with a better way to determine if you have a two
        %dimensional object
        dimrm = all(loclog,1);
        if dimrm(1)&&~emptyarr(1)
            X=ones(size(Y))*sx;
        elseif dimrm(2)&&~emptyarr(2)
            Y=ones(size(Z))*sy;
        elseif dimrm(3)&&~emptyarr(3)
            Z=ones(size(X))*sz;
        end
        [curcdata,alpha] = makecdata(V,cmap,vbound);
        hslice=surf(axh,X,Y,Z,'Cdata',curcdata,'EdgeColor','none',...
            'FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',alpha);
    elseif surftype==1
        % TODO: Make it so reshape is more generalized.
        % Making assumption that data is already in basically a flattened
        % version of the the shape of the object you want plotted.
        [curcdata,alpha] = makecdata(reshape(v,size(Xi)),cmap,vbound);
        hslice=surf(axh,Xi,Yi,Zi,'Cdata',curcdata,'EdgeColor','none',...
            'FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',alpha);
    end
else
    % reshape the data
    [X,Y,Z,V] = reshapegen(GD.dataloc,v);
    if surftype==1
        hslice=slice(axh,X,Y,Z,V,Xi,Yi,Zi);
    elseif surftype==2
        hslice=slice(axh,X,Y,Z,V,sx,sy,sz);
    end
    colormap(axh,cmap)
    caxis(vbound);
end


xlabel('\bf x [km]');
ylabel('\bf y [km]');
zlabel('\bf z [km]');
title(titlestr,'FontSize',16);
shading flat;
