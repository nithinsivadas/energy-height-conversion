function [msisData,F107A,F107,APH] = msis_irbem(timeMSIS,coords,msisVersionStr,coordType)
%% MSIS_IRBEM Produces MSIS neutral number densities and temperature in SI units. 
%The function uses IRBEM package. The input coordinates are in GDZ coordinates 
% by default (alt, lat, lon). However they can be changed to any
% that irbem allows.
% 
% Syntax:  [out,F107A,F107,APH] = msis_irbem(timeMSIS,coords)
%
% Inputs:
%    timeMSIS - An array or single element of matlab time [Nx1]
%    coords   - [X,Y,Z] coordinates [Nx3], by default GDZ coordinates
%               [alt, lat, lon] - in that order. 
%    msisVersionStr - String indicating the atmospheric model
%               'nrlmsise00','msise90','msis86'
%    coordType - Default: 0 -> GDZ (alt,lat,lon)
%               Range from 0->8 [GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH,
%               RLL], See https://svn.code.sf.net/p/irbem/code/tags/IRBEM-4.4.0/manual/user_guide.html           
%
% Outputs:
%    msisData.H - Hydrogen number density in m-3
%       -> {He,O,O2,N2,Ar,N,AnomalousOxygen} 
%                 as well in the same units
%       -> TotalMass - Total mass in kg m-3
%    F107A      - F10.7a value given as input to the MSIS model
%    F107       - F10.7 value given as input to the MSIS model
%    APH        - AP values, integrated over different time periods given
%                 as input to the MSIS model
%
% Example: 
%   [msisData,F107A,F107,APH] = msis_irbem(datenum('26 Mar 2008'),...
%       [(linspace(50,100,5))',repmat(65,5,1),repmat(-147,5,1)]);
%
% Other m-files required: f107_aph.m, onera_desp_lib_msis.m
% Subfunctions: none
% MAT-files required: none
%
% See also: msis.m, onera_desp_lib_msis.m
% Author: Nithin Sivadas, Ph.D Candidate
% Center for Space Physics, Boston University
% email address: nithin@bu.edu
% Website: 
% November 2018; Last revision: 20-Mon-2018
%------------- BEGIN CODE --------------

%%
if nargin<4
    coordType=0;
end
if nargin<3
    msisVersionStr='nrlmsise00';
end

timeaph = min(timeMSIS(:));
[F107A, F107, APH] = f107_aph(timeaph);
nCoords=size(coords,1);

% Making sure that the number of time elements, equal the number of
% coordinate enteries. 
if length(timeMSIS)==1
    timeArray = repmat(timeMSIS,nCoords,1);
else 
    timeArray = timeMSIS(:)';
    if length(timeArray) ~= nCoords
    error('The number of time elements, should match the number of coordinates');
    end
end

msisData = onera_desp_lib_msis(msisVersionStr,timeArray,coords,...
    coordType,repmat(F107A,nCoords,1),repmat(F107,nCoords,1),repmat(APH,nCoords,1));

% Converting to SI Units
msisData.He = msisData.He * 10^6; % m^-3
msisData.O = msisData.O * 10^6; % m^-3
msisData.O2 = msisData.O2 * 10^6; % m^-3
msisData.N2 = msisData.N2 * 10^6; % m^-3
msisData.Ar = msisData.Ar * 10^6; % m^-3
msisData.H = msisData.H * 10^6; % m^-3
msisData.N = msisData.N * 10^6; % m^-3
msisData.AnomalousOxygen = msisData.AnomalousOxygen * 10^6; % m^-3
msisData.TotalMass = msisData.TotalMass * 10^3; % kg/m^-3
msisData.TotalNumberDensity = msisData.He + msisData.O + msisData.N2 + ...
    msisData.Ar + msisData.H + msisData.N + msisData.O2 + ...
    msisData.AnomalousOxygen;    
msisData.Description = 'He,O,N2,O2,Ar,H,N, are all in number densities m-3; Total Mass is in Kg m-3; Temperatures in K';
end

