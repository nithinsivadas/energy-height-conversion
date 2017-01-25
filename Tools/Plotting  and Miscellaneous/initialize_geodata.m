function [] = initialize_geodata(filePath)
%initialize_geodata.m The function initializes GeoData, given the path
%where the GeoData setup.m is stored
computer=getenv('computername');
if nargin<1
    if computer=='NITHIN-SURFACE'
        filePath=['C:\Users\Nithin\Documents\GitHub\energy-height-conversion\GeoDataMATLAB\'];    
    else
        filePath=['/home/nithin/Documents/git-repos/energy-height-conversion/GeoDataMATLAB/'];
    end
end;

eval(sprintf('run %s;',([filePath,'setup.m'])));

end

