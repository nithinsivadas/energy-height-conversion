function [probes] = find_magnetic_conjunctions_at_mag_equator(time,maginput,...
    spacecraft,conjunction,stopAlt, hemiFlag, magFieldModel)
%% find_magnetic_conjunctoin.m Calculates magnetic conjunction between 
%  spacecraft, groundstation and selected conjunction probe at stopAlt heights
%--------------------------------------------------------------------------
% Input
%------
% time      : [nTx1] [u:matlabtime] An array of time points at which you need
%             to estimate magnetic conjunction
%             Default:'26 Mar 2008'
% omniData  : Data product of process_omni_data.m Contains
%             omniData.minutely.maginput, for estimate the magnetic field
%             model
% spacecraft: 
%   -> (probename).GEO: [nTx3] Time series position coordinates of a
%                       spacecraft in GEO coordinate system. There can be
%                       more than one probename. Example spacraft.thma,
%                       spacecraft.thmb and so on in the same structure. 
%   -> (probename).time: [nTx1] An array of time points corresponding to
%                        the spacecraft position data. 
%                       Example: 'pokerflat'
% conjunction:
%   -> probeStr       : [string] Select the spacecraft or ground station
%                       with respect to which the magnetic field
%                       conjunction of the other probes ought to be
%                       calculated
%                       Example: 'pokerflat' or 'thma'
%   -> radius         : [1][km] Radius of the circle within which the
%                       geodetic coordinates of the probes have to fall to
%                       be considered as magnetically conjugate with
%                       probeStr. 
% stopAlt   : [km] The altitude at which the magnetic foot point of the
%             spacecraft ought to be calculated. 
%             Default: 110km
% hemiFlag  : A flag value to specify the hemisphere of conjunction 
%                0  - same Hemisphere as start point
%               +1  - Northern Hemisphere
%               -1  - Southern Hemisphere
%                2  - opposite Hemisphere as start point
%            Default: +1
%
% magFieldModel: Flat specifying which magnetic field to use. Default:11
%                 0   = no external field 
%                 1   = Mead & Fairfield [1975] 
%                     (uses 0?Kp?9 - Valid for rGEO?17. Re) 
%                 2   = Tsyganenko short [1987] 
%                     (uses 0?Kp?9 - Valid for rGEO?30. Re) 
%                 3   = Tsyganenko long [1987] 
%                     (uses 0?Kp?9 - Valid for rGEO?70. Re) 
%                 4   = Tsyganenko [1989c] 
%                     (uses 0?Kp?9 - Valid for rGEO?70. Re) 
%                 5   = Olson & Pfitzer quiet [1977] 
%                     (default - Valid for rGEO?15. Re) 
%                 6   = Olson & Pfitzer dynamic [1988] 
%                     (uses 5.?dens?50., 300.?velo?500., 
%                     -100.?Dst?20. - Valid for rGEO?60. Re) 
%                 7   = Tsyganenko [1996] 
%                     (uses -100.?Dst (nT)?20., 0.5?Pdyn (nPa)?10., 
%                     |ByIMF| (nT)?10., |BzIMF| (nT)?10. - Valid for rGEO?40. Re) 
%                 8   = Ostapenko & Maltsev [1997] 
%                     (uses dst,Pdyn,BzIMF, Kp) 
%                 9   = Tsyganenko [2001] 
%                     (uses -50.?Dst (nT)?20., 0.5?Pdyn (nPa)?5., |ByIMF| (nT)?5., 
%                     |BzIMF| (nT)?5., 0.?G1?10., 0.?G2?10. - Valid for xGSM?-15. Re) 
%                 10 =Tsyganenko [2001] storm  
%                     (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 
%                     - there is no upper or lower limit 
%                     for those inputs - Valid for xGSM?-15. Re) 
%                 11 =Tsyganenko [2004] storm  
%                     (uses Dst, Pdyn, ByIMF, BzIMF, 
%                     W1, W2, W3, W4, W5, W6 
%                     - there is no upper or lower limit for those inputs 
%                     - Valid for xGSM?-15. Re) 
%                 12 =Alexeev [2000] 
%                     - also known as Paraboloid model - 
%                     Submitted to ISO  (uses Dens, velo, Dst, BzIMF, AL) 
%
%--------------------------------------------------------------------------
% Output
%-------
% probes:
% -> flag : [nTxnProbes] Each element in this array corresponds to a time
%        point and probe. If the value is +1, at that time point that
%        particular probe is in magnetic conjunction as defined by the
%        conjunction.radius with conjunction.probeStr. If it is 0, then
%        there is no conjunction.
%   -> probeNames: [1xnProbes][cell] Each element corresponds to column in flag, 
%                  and is the name of the probe that it corresponds to.  
%   -> GDZ       : [nT x 2] Lat and Lon coordinates of the probe for
%                  the correspoinding input times
%   -> GEO       : [nT x 3] X,Y,Z- GEO cartesian coordinates of probes 
%                  which are spacecrafts
%--------------------------------------------------------------------------
% Modified: 13th Feb 2017 
% Created : 13th Feb 2017
% Author  : Nithin Sivadas
% Ref     : 
% Notes   : The CDF file is stored in an awkard location - under each event
%           date. Only processed .mat files need to be stored here. Change
%           this.
%--------------------------------------------------------------------------
%% IRBEM options
options = [0,0,0,0,0];
sysaxes = 1; %GEO Input coordinates
if nargin<8
    magFieldModel = 7; % External magnetic field model: TSY 1996
end
if nargin<7
    hemiFlag = +1; % Northern hemisphere conjunctions
end
if nargin<6
    stopAlt = 110; %Visible auroral altitude
end
%% Calculating magnetic foot prints
nTime = length(time);
probeNames = fieldnames(spacecraft);
nProbes = length(probeNames);
multiWaitbar('find_magnetic_conjunctions',0);
for thisSC = 1:1:nProbes
    thisProbe = char(probeNames(thisSC));
    multiWaitbar('find_magnetic_conjunctions','Increment',1/nProbes,'Relabel',['Magnetic conjunctions for ',thisProbe]);
    GEO=interp1(spacecraft.(thisProbe).time,spacecraft.(thisProbe).GEO,time,'linear');
    GDZNorth=onera_desp_lib_find_foot_point(magFieldModel,options,...
        sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),stopAlt,+1,maginput); % alt lat lon
%     GDZSouth=onera_desp_lib_find_foot_point(magFieldModel,options,...
%         sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),stopAlt,-1,maginput); % alt lat lon
%     [Bmin,magEqPointGEO]=onera_desp_lib_find_magequator(magFieldModel,options,...
%         sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),maginput); % alt lat lon
    [BGEO]=onera_desp_lib_get_field(magFieldModel,options,...
        sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),maginput); % alt lat lon
    
%     [Lm,Lstar,Blocal,Bmin,XJ,MLT] = onera_desp_lib_make_lstar(magFieldModel,options,...
%         sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),maginput);
%     [xHEMI] = onera_desp_lib_get_hemi(magFieldModel,options,sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),maginput);
    %     [Lm,Blocal,Bmin,J,POSIT]=onera_desp_lib_trace_field_line_towards_earth...
%         (magFieldModel,options,sysaxes,time,GEO(:,1),GEO(:,2),GEO(:,3),maginput);
    multiWaitbar('find_foot_point','Reset');
%     disp('end IRBEM\n');
    probes.(thisProbe).GDZ = GDZNorth(:,2:3); %lat lon
%     probes.(thisProbe).GDZNorth = GDZNorth(:,2:3); %lat lon
%     probes.(thisProbe).GDZSouth = GDZSouth(:,2:3); %lat lon
    probes.(thisProbe).GEO = GEO; %lat lon
%     probes.(thisProbe).magEqGEO = magEqPointGEO;
%     probes.(thisProbe).Bmin = Bmin; % minimum value of B-field
%     probes.(thisProbe).Lm = Lm;
%     probes.(thisProbe).Lstar = Lstar;
    probes.(thisProbe).BlocalGEO = BGEO;
%     probes.(thisProbe).MLT = MLT;
%     probes.(thisProbe).xHEMI = xHEMI; %+1 - North mag. hemisphere, -1 South mag. hemisphere, 0 cannot be defined
    multiWaitbar(['Magnetic conjunctions for ',thisProbe],'Relabel','find_magnetic_conjunctions');
end
multiWaitbar('find_magnetic_conjunctions','Reset');
%% Identifying conjunction
% GDZ=repmat(groundstation.GDZ,nTime,1);
% probes.(groundstation.field).GDZ = GDZ; %lat lon
fieldNames = fieldnames(probes);
indx = find(strcmp(fieldNames,conjunction.probeStr));
conjunctionCoords.GDZ = probes.(char(fieldNames(indx))).GDZ;
conjunctionCoords.radius = repmat(conjunction.radius,nTime,1);
flag = find_geodetic_conjunction(probes,conjunctionCoords,stopAlt); %   can add more decision parameters in future
% probeNames(thisSC+1)={groundstation.field};
probes.probeNames = probeNames';
probes.flag = flag;
end
