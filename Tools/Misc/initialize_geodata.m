function [] = initialize_geodata(filePath)
%initialize_geodata.m The function initializes GeoData, given the path
%where the GeoData setup.m is stored
rootPathStr=initialize_root_path;
eval(sprintf('run %s;',([rootPathStr,'energy-height-conversion\GeoDataMATLAB\setup.m'])));
end

