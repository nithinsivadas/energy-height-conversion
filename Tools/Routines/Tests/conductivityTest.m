alt = 60:1:500;
lat = 65.1;
lon = -147.5;
time = datenum('24 Jan 2012 11:19');

get_conductivity_v2(alt,[],lat,lon,time,1,{'all'},true);

%%
% data = iri2016f90(time,alt,lat,lon,true);