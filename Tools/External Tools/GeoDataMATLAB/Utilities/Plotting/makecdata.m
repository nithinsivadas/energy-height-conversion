function varargout = makecdata(object,cmap,bounds)
% by John Swoboda
% This function will take an array of any type and turn it into to colormap
% data for an image. The final array will have the same dimensional
% structure except the last dimension will be the colormap data.
% Inputs
% object - This is an n-dimensional array that will be plotted. 
% cmap - This is a colormap that will be applied to the data. This will be
% a ncolorx3 array.
% bounds - This is the bounds of color scale of the array. 
% Outputs
% Cdata- The final color data. This is the same shape as object but with
% one additional dimension.

% Flattent the array
oldshape = size(object);
object= object(:);
oldobj = object;
% Find the minimum and span of the colors.
minb = bounds(1);
maxb = bounds(2);
cspan = maxb-minb;
[ncolors,ndim] = size(cmap);
% saturate bottom and top
object = min(object,maxb);
object = max(object,minb);
object(isnan(oldobj)) = minb;
object = ceil((ncolors-1)*(object-minb)/cspan)+1;

alphamap = ~isnan(oldobj);
alphamap = double(reshape(alphamap,oldshape));
Cdata = cmap(object,:);
Cdata = reshape(Cdata,[oldshape,ndim]);

varargout{1} = Cdata;

if nargout >1
    varargout{2} = alphamap;
end
    