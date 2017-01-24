function [] = include_new_colormaps(filePath)
%include_new_colormaps.m The function initializes a few colormaps, given the path
%where the colormaps are stored

if nargin<1
    filePath=['/home/nithin/Documents/git-repos/energy-height-conversion/Tools/External Tools/Colormaps'];
end;

eval(sprintf('addpath(''%s'');',([filePath])));

end
