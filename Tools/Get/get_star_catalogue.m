function [data] = get_star_catalogue(fitsFile)
%get_star_catalogue Reads the hipparcos star catalogue, 
%                   Give ICRS RA, and DEC in degrees. 
%   Detailed explanation goes here
% Use this: http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=V/137D
% to get Hipparcus data
%
%  Find with star name:
%  starIndx = find(contains(deblank(data.name),'Vega','IgnoreCase',true));
%  deblank(data.name(starIndx))
%
if nargin<1
    fitsFile = [initialize_root_path,filesep,'energy-height-conversion',...
        filesep,'Tools',filesep,'External Tools',filesep,...
       'skymap',filesep,'hipparcos_extended_catalogue_J2000.fit'];
end

temp=fitsread(fitsFile,'binarytable',1);
temp2=fitsread(fitsFile,'binarytable',2);

data.HIP = temp{3};
data.vmag = temp{10};
data.luminosity = temp{11};
data.RA = deg2rad(temp{1});
data.DEC = deg2rad(temp{2});
data.parallax = temp{6}./1000;
data.pmDEC = deg2rad(temp{8}./3600000);
data.pmRA = deg2rad(temp{7}./3600000)./cosd(temp{5});
data.RV = temp{9};
data.name = temp{12};
data.angularRadius = temp2{8};

m1 = min(data.vmag);
data.relIntensity = 10.^(0.4.*(m1-data.vmag));

data.description.HIP = {'Identifier (HIP number) (H1)','Unit: #'};
data.description.vmag = {'Magnitude in Johnson V (H5)','Unit: magnitude'};
data.description.luminosity = {'Stellar luminosity','Unit: Lsun solar luminosity'};
data.description.RA = {'Right ascension (KF5) at Epoch=J2000','Unit: rad'};
data.description.DEC = {'Declination (FK5) at Epoch=J2000','Unit: rad'};
data.description.parallax = {'Trigonometric parallax','Unit: arcsec'};
data.description.pmRA = {'Proper motion in RA','Unit: rad per year (rad/year)'};
data.description.pmDEC = {'Proper motion in Declination','Unit: rad per year (rad/year)'};
data.description.RV = {'Radial Velocity','Unit: km/s'};
data.description.name = {'Star Name'};
data.description.angularRadius = {'Angular Radius 1sigma','Unit: deg'};
data.description.relIntensity = {'Relative intenisty of stars','I2/I1 - calculated from 10^0.4(m1-m2)'};


end

