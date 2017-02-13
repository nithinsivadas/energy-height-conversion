function [] = include_new_colormaps(filePath)
%% include_new_colormaps.m The function initializes a few colormaps, given the path of the color maps
rootPathStr = initialize_root_path();
if nargin<1
    filePath=[rootPathStr,'/Tools/External Tools/Colormaps'];
end;

eval(sprintf('addpath(''%s'');',([filePath])));

end
