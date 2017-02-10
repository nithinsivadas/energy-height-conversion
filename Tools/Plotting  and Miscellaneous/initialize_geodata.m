function [] = initialize_geodata(filePath)
%initialize_geodata.m The function initializes GeoData, given the path
%where the GeoData setup.m is stored
rootPathStr=root_path;
eval(sprintf('run %s;',([rootPathStr,'GeoDataMATLAB\setup.m'])));
end

