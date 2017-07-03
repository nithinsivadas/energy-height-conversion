function [ data ] = get_all_time_series_data()
% get_time_series_data Provides a structure Data, that contains all the 
% data files that can be later be used for plotting.

% In this version of the function, you get the following data in the
% structures:

% Output:
% data(1): IMF Bz [1-D Time Series, with 1 variable]
% data(2): AL & AU indices [1-D Time Series, with 2 variables]
% data(3): ASI Keogram of FYKN [2-D Time Series]
% data(4): Electron density from PFISR [2-D Time Series]
% data(5): Ionosphere Inverse Energy Spectra (Vickrey Model) [2-D Time Series]
% data(6): Ionosphere Inverse Energy Spectra (SIC Model) [2-D Time Series]
% data(7): B field Z and X from FGM [1-D Time Series, with 2 variables]
% data(8): Energy spectra from thd [2-D Time Series]
% data(9): Ion Anisotropy [1-D Time Series, with 2 variables]
% data(10): Electron Anisotropy [1-D Time Series, with 2 variables]
% data(11): Wave power 1-1000Hz from FBK [2-D Time Series]
% data(12): Wave power <1Hz from SCM along X [2-D Time Series]
% data(13): Wave power <1Hz from EFI along Y [2-D Time Series]
% data(14): Wave power <1Hz from FGM along X [2-D Time Series]
% data(15): ASI Keogram of GAKO [2-D Time Series]
% data(16): Cumulative Inverse Energy Spectra (Vickrey Model) [2-D Time Series]
% data(17): Cumulative Inverse Energy Spectra (SIC Model) [2-D Time Series]
% data(18): Cumulative Energy spectra from thd [2-D Time Series]
% data(19): Parallel Potential [1-D Time Series, with 1 Variable]
%----------------------------------------------------------------------------
% Modified: 21st Sep 2016| 25th Jan 2017 
% Created : 21st Sep 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------

%% Loading variables

computer=getenv('computername');
if computer=='NITHIN-SURFACE'
load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Ground\pfisr_data_26_Mar_2008.mat')
load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Ground\thmasi_data_26_Mar_2008.mat')

load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Space\thd_data_26_Mar_2008.mat')
load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Space\parallel_potential.mat')

load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Global\IMF_data_26_Mar_2008.mat')
load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Global\GeomagneticIndex_data_26_Mar_2008.mat')    
else
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Ground/pfisr_data_26_Mar_2008.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Ground/thmasi_data_26_Mar_2008.mat')

load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/parallel_potential.mat')

load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Global/IMF_data_26_Mar_2008.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Global/GeomagneticIndex_data_26_Mar_2008.mat')
end
%% Organizing the data accoring to number

% Global Measurements
% %  IMF  Bz
% thisData=1;
% Data(thisData).time = imf_ace.time;
% Data(thisData).yAxis = imf_ace.Bz;
% Data(thisData).label = ({'IMF B_z','[nT]'});
% Data(thisData).y.tick=[-10, -5, 0, 5, 10];
% Data(thisData).y.tickLabel={'-10', '-5', '0', '5', '10'};
% Data(thisData).y.lim = [min(Data(thisData).yAxis)-5 max(Data(thisData).yAxis)];
% Data(thisData).legend = ['B_z'];
%  IMF_OMNI  Bz
thisData=1;
[data(thisData).yAxis, data(thisData).time] = crop_time(imf_omni.Bz', imf_omni.time',...
    datenum('26 Mar 2008 08:00'), datenum('26 Mar 2008 13:00'));
data(thisData).label = ({'IMF B_z','[nT]'});
data(thisData).y.tick=[-10, -5, 0, 5, 10];
data(thisData).y.tickLabel={'-10', '-5', '0', '5', '10'};
data(thisData).y.lim = [min(data(thisData).yAxis)-5 max(data(thisData).yAxis)];
data(thisData).legend = ['B_z'];
% Geomagnetic Index (AL & AU)
thisData=2;
data(thisData).time = geomagIndex.time';
data(thisData).yAxis = geomagIndex.al;
data(thisData).yAxis2 = geomagIndex.au;
data(thisData).label = ({'Mag Activity','[\muT]'});
data(thisData).y.tick=[-1500, -1000, -500, 0, 500];
data(thisData).y.tickLabel={'-1.5', '-1.0', '-0.5', '0', '0.5'};
data(thisData).y.lim = [min(data(thisData).yAxis)-200 max(data(thisData).yAxis2+200)];
data(thisData).legend = [{'AL'}, {'AU'}];
% Ground Measurements
% FYKN ASI Keogram 
thisData=3;
data(thisData).time = asi.fykn.time;
data(thisData).yAxis = asi.fykn.mlat;
data(thisData).zValue = asi.fykn.intensity';
data(thisData).label = ({'FYKN 67.22^0N','Keogram','[MLAT]'});
data(thisData).y.tick=[60, 62, 64, 66, 68, 70];
data(thisData).y.tickLabel={'60^0N', '62^0N', '64^0N', '66^0N', '68^0N', '70^0N'};
data(thisData).color.label={'Intensity'};
data(thisData).y.lim=[60 70];
data(thisData).color.lim=10.^[3.4 3.8];
data(thisData).color.tick=10.^[3.4 3.5 3.6 3.7 3.8];
data(thisData).color.tickLabel={'10^3^.^4','10^3^.^5', '10^3^.^6', '10^3^.^7', '10^3^.^8'};


% Electron density from PFISR 
thisData=4;
data(thisData).time = pfisrBeamAvg.density.time;
data(thisData).yAxis = pfisrBeamAvg.density.alt;
data(thisData).zValue = pfisrBeamAvg.density.electron;
data(thisData).label = ({'PFISR', 'N_e Density','Altitude [km]'});
data(thisData).color.lim=10.^[8 12];
data(thisData).color.tick=10.^[8 9 10 11 12];
data(thisData).color.tickLabel={'10^8','10^9', '10^1^0', '10^1^1', '10^1^2'};
data(thisData).color.label={'[m^-^3]'};
data(thisData).y.tick=[70, 80, 90, 100, 120, 140];
data(thisData).y.tickLabel={'70', '80', '90', '100', '120', '140'};
data(thisData).y.lim=[60 150];

% Ionosphere Inverse Energy Spectra (Vickrey Model)
thisData=5;
data(thisData).time = pfisrBeamAvg.energySpectraVickery.time;
data(thisData).yAxis = pfisrBeamAvg.energySpectraVickery.energyBin;
data(thisData).zValue = pfisrBeamAvg.energySpectraVickery.energyFlux;
data(thisData).label = ({'PFISR Inversion','(Vickrey et.al 1982)','Diff. energy flux','[keV]'});
data(thisData).color.lim=10.^[6 12];
data(thisData).color.tick=10.^[6 8 10 12];
data(thisData).color.tickLabel={'10^6','10^8', '10^1^0', '10^1^2'};
data(thisData).color.label={'[eV/m^2 sr s eV]'};
data(thisData).y.tick=10.^[3 4 5 5.477];
data(thisData).y.tickLabel={'1','10','100','300'};
data(thisData).y.lim=[1*10^3 900*10^3];

% Ionosphere Inverse Energy Spectra (SIC Model)
thisData=6;
data(thisData).time = pfisrBeamAvg.energySpectraSIC.time';
data(thisData).yAxis = pfisrBeamAvg.energySpectraSIC.energyBin;
data(thisData).zValue = pfisrBeamAvg.energySpectraSIC.energyFlux;
data(thisData).label = ({'PFISR Inversion','SIC Model','Diff. energy flux','[keV]'});
data(thisData).color.lim=10.^[6 12];
data(thisData).color.tick=10.^[6 8 10 12];
data(thisData).color.tickLabel={'10^6','10^8', '10^1^0', '10^1^2'};
data(thisData).color.label={'[eV/m^2 sr s eV]'};
data(thisData).y.tick=10.^[3 4 5 5.477];
data(thisData).y.tickLabel={'1','10','100','300'};
data(thisData).y.lim=[1*10^3 900*10^3];

% Space Measurements
% B field Z and X from FGM
thisData=7;
data(thisData).time = thd.magneticVectorFieldFGM.time;
data(thisData).yAxis = thd.magneticVectorFieldFGM.Bz;
data(thisData).yAxis2 = thd.magneticVectorFieldFGM.Bx;
data(thisData).label = ({'Thm-D','B-field','[nT]'});
data(thisData).y.lim = [min(data(thisData).yAxis2)-10 max(data(thisData).yAxis)+5];
data(thisData).y.tick=[-75, -50, -25, 0, 25, 50, 75];
data(thisData).y.tickLabel={'-75', '-50', '-25', '0', '25', '50', '75'};
data(thisData).legend = [{'B_z'},{'B_x'}];

% Energy spectra from thd
thisData=8;
data(thisData).time = thd.energySpectra.time;
data(thisData).yAxis = thd.energySpectra.energyBin;
data(thisData).zValue = thd.energySpectra.energyFlux';
data(thisData).label = ({'Thm-D','Diff. e^- energy flux','[keV]'});
data(thisData).color.lim=10.^[6 12];
data(thisData).color.tick=10.^[6 8 10 12];
data(thisData).color.tickLabel={'10^6','10^8', '10^1^0', '10^1^2'};
data(thisData).color.label={'[eV/m^2 sr s eV]'};
data(thisData).y.tick=10.^[3 4 5 5.477];
data(thisData).y.tickLabel={'1','10','100','300'};
data(thisData).y.lim=[1*10^3 900*10^3];

% Ion Anisotropy
thisData=9;
data(thisData).time = thd.pitchAngleAnisotropy.sst.ionTime;
data(thisData).yAxis = thd.pitchAngleAnisotropy.sst.ionAnisotropy;
data(thisData).time2 = thd.pitchAngleAnisotropy.esa.ionTime;
data(thisData).yAxis2 = thd.pitchAngleAnisotropy.esa.ionAnisotropy;
data(thisData).label = ({'Thm-D',' i^+ Anisotropy','[a.u.]'});
data(thisData).y.tick=[-1.0 -0.5 0 0.5 1.0];
data(thisData).y.tickLabel={'-1.0','-0.5','0','0.5' , '1.0'};
data(thisData).y.lim = [-1.0 +1.0];
data(thisData).legend = [{'> 25 keV'}, {'< 25 keV'}];
% Electron Anisotropy
thisData=10;
data(thisData).time = thd.pitchAngleAnisotropy.sst.electronTime;
data(thisData).yAxis = thd.pitchAngleAnisotropy.sst.electronAnisotropy;
data(thisData).time2 = thd.pitchAngleAnisotropy.esa.electronTime;
data(thisData).yAxis2 = thd.pitchAngleAnisotropy.esa.electronAnisotropy;
data(thisData).label = ({'Thm-D','e^- Anisotropy','[a.u.]'});
data(thisData).y.lim = [-1.0 +1.0];
data(thisData).y.tick=[-1.0, -0.5, 0, 0.5, 1.0];
data(thisData).y.tickLabel={'-1.0','-0.5','0','0.5' , '1.0'};
data(thisData).legend = [{'> 30 keV'}, {'< 30 keV'}];

% Wave power 1-1000Hz from FBK
thisData=11;
data(thisData).time = thd.waveHighFrequencyFB_EFI.time;
data(thisData).yAxis = thd.waveHighFrequencyFB_EFI.frequency;
data(thisData).zValue = thd.waveHighFrequencyFB_EFI.power';

% Correcting for 0 values:
minzValue=min(min(((data(thisData).zValue(data(thisData).zValue>0)))));
data(thisData).zValue(data(thisData).zValue==0)=minzValue;

data(thisData).label = ({'Thm-D','FB','[Hz]'});
data(thisData).color.label={'[mV m^-^1 Hz^-^0^.^5]'};
data(thisData).y.tick=10.^[0 1 2 3];
data(thisData).y.tickLabel={'1','10','100','10^3'};
data(thisData).y.lim=[min(data(thisData).yAxis) max(data(thisData).yAxis)];
data(thisData).color.lim=10.^[-2 1];
data(thisData).color.tick=10.^[-2 -1 0 1];
data(thisData).color.tickLabel={'10^-^2','10^-^1', '10^0', '10^1'};

% Wave power <1Hz from SCM along X
thisData=12;
data(thisData).time = thd.waveLowFrequencySCM.time;
data(thisData).yAxis = thd.waveLowFrequencySCM.frequency;
data(thisData).zValue = thd.waveLowFrequencySCM.powerGsmX';

% Correcting for 0 values:
minzValue=min(min(((data(thisData).zValue(data(thisData).zValue>0)))));
data(thisData).zValue(data(thisData).zValue==0)=minzValue;

data(thisData).label = ({'Thm-D','SCM GSM-X','[Hz]'});
data(thisData).color.label={'[nT Hz^-^0^.^5]'};
data(thisData).y.tick=10.^[-2 -1 0];
data(thisData).y.lim=[min(data(thisData).yAxis) max(data(thisData).yAxis)];
data(thisData).y.tickLabel={'10^-^2','10^-^1','10^0'};
data(thisData).color.lim=10.^[-4 3];
data(thisData).color.tick=10.^[-4 -2 0 2];
data(thisData).color.tickLabel={'10^-^4','10^-^2', '10^0', '10^2'};

% Wave power <1Hz from EFI along Y
thisData=13;
data(thisData).time = thd.waveLowFrequencyEFI.time;
data(thisData).yAxis = thd.waveLowFrequencyEFI.frequency;
data(thisData).zValue = thd.waveLowFrequencyEFI.wavePowerGsmY';

% Correcting for 0 values:
minzValue=min(min(((data(thisData).zValue(data(thisData).zValue>0)))));
data(thisData).zValue(data(thisData).zValue==0)=minzValue;

% Calculating Ion Gyro Freqyencies
B = interp1(thd.magneticVectorFieldFGM.time,thd.magneticVectorFieldFGM.BTotal,data(thisData).time).*10^-9; % [T]
e = 1.602*10^-19; % [Coulombs]
m = 1.672619*10^-27; % [kg] 
f_HI = e*B/(2*pi*m); % Hz
f_OI = e*B/(2*pi*16*m); % Hz
f_HeI = e*B/(2*pi*4*m); % Hz
f_HeII = 2*e*B/(2*pi*4*m); %Hz
data(thisData).label = ({'Thm-D','EFI GSM-Y','[Hz]'});
data(thisData).color.label={'[mV m^-^1 Hz^-^0^.^5]'};
data(thisData).y.tick=10.^[-2 -1 0];
data(thisData).y.tickLabel={'10^-^2','10^-^1','10^0'};
data(thisData).y.lim=[min(data(thisData).yAxis) max(data(thisData).yAxis)];
data(thisData).color.lim=10.^[-4 3];
data(thisData).color.tick=10.^[-4 -2 0 2];
data(thisData).color.tickLabel={'10^-^4','10^-^2', '10^0', '10^2'};
data(thisData).y.f_HI =f_HI;
data(thisData).y.f_OI =f_OI;
data(thisData).y.f_HeI =f_HeI;
data(thisData).y.f_HeII =f_HeII;

% Wave power <1Hz from FGM along X
thisData=14;
data(thisData).time = thd.waveLowFrequencyFGM.time;
data(thisData).yAxis = thd.waveLowFrequencyFGM.frequency;
data(thisData).zValue = thd.waveLowFrequencyFGM.powerGsmX';

% Correcting for 0 values:
minzValue=min(min(((data(thisData).zValue(data(thisData).zValue>0)))));
data(thisData).zValue(data(thisData).zValue==0)=minzValue;

data(thisData).label = ({'Thm-D','FGM GSM-X','[Hz]'});
data(thisData).color.label={'[nT Hz^-^0^.^5]'};
data(thisData).y.tick=10.^[-2 -1 0];
data(thisData).y.tickLabel={'10^-^2','10^-^1','10^0'};
data(thisData).y.lim=[min(data(thisData).yAxis) max(data(thisData).yAxis)];
data(thisData).color.lim=10.^[-4 3];
data(thisData).color.tick=10.^[-4 -2 0 2];
data(thisData).color.tickLabel={'10^-^4','10^-^2', '10^0', '10^2'};

% Additional Ground Measurements
% GAKO ASI Keogram 
thisData=15;
data(thisData).time = asi.gako.time;
data(thisData).yAxis = asi.gako.mlat;
data(thisData).zValue = asi.gako.intensity';
data(thisData).label = ({'GAKO 63.05^0N','Keogram','[MLAT]'});
data(thisData).y.tick=[60, 62, 64, 66, 68, 70];
data(thisData).y.tickLabel={'60^0N', '62^0N', '64^0N', '66^0N', '68^0N', '70^0N'};
data(thisData).color.label={'Intensity'};
data(thisData).y.lim=[60 70];
data(thisData).color.lim=10.^[3.4 3.8];
data(thisData).color.tick=10.^[3.4 3.5 3.6 3.7 3.8];
data(thisData).color.tickLabel={'10^3^.^4','10^3^.^5', '10^3^.^6', '10^3^.^7', '10^3^.^8'};


% Cumulative Energy Fluxes
% PFISR Vickery
thisData=16;
data(thisData).zValue= diff_to_cumu_flux(...,
    pfisrBeamAvg.energySpectraVickery.energyFlux,...
    pfisrBeamAvg.energySpectraVickery.energyBin);
data(thisData).zValue=100*data(thisData).zValue/max(max(data(thisData).zValue));
data(thisData).time = pfisrBeamAvg.energySpectraVickery.time;
data(thisData).yAxis = pfisrBeamAvg.energySpectraVickery.energyBin;

halfEnergyIndex = find_median_energy( data(thisData).zValue, 1 );
data(thisData).yValue = data(thisData).yAxis(halfEnergyIndex)'; % 50% energy flux Energy

data(thisData).label = ({'PFISR Inversion','(Vickrey et.al 1982)','Norm. cumulative','energy flux','[keV]'});
data(thisData).color.lim=10.^[0 2];
data(thisData).color.tick=10.^[0 1 2];
data(thisData).color.tickLabel={'0','10', '100'};
data(thisData).color.label={'[%]'};
data(thisData).y.tick=10.^[3 4 5 5.477];
data(thisData).y.tickLabel={'1','10','100','300'};
data(thisData).y.lim=[1*10^3 900*10^3];


% PFISR (SIC Model)
thisData=17;
data(thisData).time = pfisrBeamAvg.energySpectraSIC.time';
data(thisData).yAxis = pfisrBeamAvg.energySpectraSIC.energyBin;
data(thisData).zValue= diff_to_cumu_flux(...,
    pfisrBeamAvg.energySpectraSIC.energyFlux,...
    pfisrBeamAvg.energySpectraSIC.energyBin);
data(thisData).zValue=100*data(thisData).zValue/max(max(data(thisData).zValue));

halfEnergyIndex = find_median_energy( data(thisData).zValue, 1 );
data(thisData).yValue = data(thisData).yAxis(halfEnergyIndex)'; % 50% energy flux Energy

data(thisData).label = ({'PFISR Inversion','SIC Model','Norm. cumulative','energy flux','[keV]'});
data(thisData).color.lim=10.^[0 2];
data(thisData).color.tick=10.^[0 1 2];
data(thisData).color.tickLabel={'0','10', '100'};
data(thisData).color.label={'[%]'};
data(thisData).y.tick=10.^[3 4 5 5.477];
data(thisData).y.tickLabel={'1','10','100','300'};
data(thisData).y.lim=[1*10^3 900*10^3];


% THEMIS
thisData=18;
data(thisData).time = thd.energySpectra.time;
data(thisData).yAxis = thd.energySpectra.energyBin;
data(thisData).zValue= diff_to_cumu_flux(...,
    thd.energySpectra.energyFlux',...
    thd.energySpectra.energyBin');
data(thisData).zValue=100*data(thisData).zValue/max(max(data(thisData).zValue));

halfEnergyIndex = find_median_energy( data(thisData).zValue, 1 );
data(thisData).yValue = data(thisData).yAxis(halfEnergyIndex)'; % 50% energy flux Energy

data(thisData).label = ({'Thm-D','Norm. cumulative','energy flux','[keV]'});
data(thisData).color.lim=10.^[0 2];
data(thisData).color.tick=10.^[0 1 2];
data(thisData).color.tickLabel={'0','10', '100'};
data(thisData).color.label={'[%]'};
data(thisData).y.tick=10.^[3 4 5 5.477];
data(thisData).y.tickLabel={'1','10','100','300'};
data(thisData).y.lim=[1*10^3 900*10^3];

% Parallel Potential
% by estimating loss cone angle flux @ ~1.5deg in plasma sheet 
thisData=19;
data(thisData).time = parallel_potential.time;
data(thisData).yAxis = parallel_potential.V;
temp=interp1(data(16).time, data(16).yValue, data(18).time);
data(thisData).yAxis2 = (temp-data(18).yValue)/1000;
data(thisData).label = ({'\Delta \phi_|_|','[kV]'});
data(thisData).y.tick=[-20, 0, 20, 40, 80];
data(thisData).y.tickLabel={'-20','0', '20','40','80'};
data(thisData).y.lim = [-20 +80];
data(thisData).legend = [{'Difference of Peaks'},{'Difference of Median Energy'}];


end

