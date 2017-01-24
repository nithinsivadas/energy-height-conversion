function [] = initialize_geodata(filePath)
%initialize_geodata.m The function initializes GeoData, given the path
%where the GeoData setup.m is stored

if nargin<1
    filePath=['/home/nithin/Documents/git-repos/energy-height-conversion/GeoDataMATLAB/'];
end;

eval(sprintf('run %s;',([filePath,'setup.m'])));

end

