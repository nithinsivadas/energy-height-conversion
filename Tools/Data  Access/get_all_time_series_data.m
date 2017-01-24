function [ Data ] = get_all_time_series_data()
% get_time_series_data Provides a structure Data, that contains all the 
% data files that can be later be used for plotting.

% In this version of the function, you get the following data in the
% structures:

% Output:
% Data(1): IMF Bz [1-D Time Series, with 1 variable]
% Data(2): AL & AU indices [1-D Time Series, with 2 variables]
% Data(3): ASI Keogram of FYKN [2-D Time Series]
% Data(4): Electron density from PFISR [2-D Time Series]
% Data(5): Ionosphere Inverse Energy Spectra (Vickrey Model) [2-D Time Series]
% Data(6): Ionosphere Inverse Energy Spectra (SIC Model) [2-D Time Series]
% Data(7): B field Z and X from FGM [1-D Time Series, with 2 variables]
% Data(8): Energy spectra from thd [2-D Time Series]
% Data(9): Ion Anisotropy [1-D Time Series, with 2 variables]
% Data(10): Electron Anisotropy [1-D Time Series, with 2 variables]
% Data(11): Wave power 1-1000Hz from FBK [2-D Time Series]
% Data(12): Wave power <1Hz from SCM along X [2-D Time Series]
% Data(13): Wave power <1Hz from EFI along Y [2-D Time Series]
% Data(14): Wave power <1Hz from FGM along X [2-D Time Series]
% Data(15): ASI Keogram of GAKO [2-D Time Series]
% Data(16): Cumulative Inverse Energy Spectra (Vickrey Model) [2-D Time Series]
% Data(17): Cumulative Inverse Energy Spectra (SIC Model) [2-D Time Series]
% Data(18): Cumulative Energy spectra from thd [2-D Time Series]
% Data(19): Parallel Potential [1-D Time Series, with 1 Variable]

%% Loading variables
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Ground/pfisr_data_26_Mar_2008.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Ground/thmasi_data_26_Mar_2008.mat')

load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/parallel_potential.mat')

load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Global/IMF_data_26_Mar_2008.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Global/GeomagneticIndex_data_26_Mar_2008.mat')

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
[Data(thisData).yAxis, Data(thisData).time] = time_crop(imf_omni.Bz', imf_omni.time',...
    datenum('26 Mar 2008 08:00'), datenum('26 Mar 2008 13:00'));
Data(thisData).label = ({'IMF B_z','[nT]'});
Data(thisData).y.tick=[-10, -5, 0, 5, 10];
Data(thisData).y.tickLabel={'-10', '-5', '0', '5', '10'};
Data(thisData).y.lim = [min(Data(thisData).yAxis)-5 max(Data(thisData).yAxis)];
Data(thisData).legend = ['B_z'];
% Geomagnetic Index (AL & AU)
thisData=2;
Data(thisData).time = geomagIndex.time';
Data(thisData).yAxis = geomagIndex.al;
Data(thisData).yAxis2 = geomagIndex.au;
Data(thisData).label = ({'Mag Activity','[\muT]'});
Data(thisData).y.tick=[-1500, -1000, -500, 0, 500];
Data(thisData).y.tickLabel={'-1.5', '-1.0', '-0.5', '0', '0.5'};
Data(thisData).y.lim = [min(Data(thisData).yAxis)-200 max(Data(thisData).yAxis2+200)];
Data(thisData).legend = [{'AL'}, {'AU'}];
% Ground Measurements
% FYKN ASI Keogram 
thisData=3;
Data(thisData).time = asi.fykn.time;
Data(thisData).yAxis = asi.fykn.mlat;
Data(thisData).zValue = asi.fykn.intensity';
Data(thisData).label = ({'FYKN 67.22^0N','Keogram','[MLAT]'});
Data(thisData).y.tick=[60, 62, 64, 66, 68, 70];
Data(thisData).y.tickLabel={'60^0N', '62^0N', '64^0N', '66^0N', '68^0N', '70^0N'};
Data(thisData).color.label={'Intensity'};
Data(thisData).y.lim=[60 70];
Data(thisData).color.lim=10.^[3.4 3.8];
Data(thisData).color.tick=10.^[3.4 3.5 3.6 3.7 3.8];
Data(thisData).color.tickLabel={'10^3^.^4','10^3^.^5', '10^3^.^6', '10^3^.^7', '10^3^.^8'};


% Electron density from PFISR 
thisData=4;
Data(thisData).time = pfisrBeamAvg.density.time;
Data(thisData).yAxis = pfisrBeamAvg.density.alt;
Data(thisData).zValue = pfisrBeamAvg.density.electron;
Data(thisData).label = ({'PFISR', 'N_e Density','Altitude [km]'});
Data(thisData).color.lim=10.^[8 12];
Data(thisData).color.tick=10.^[8 9 10 11 12];
Data(thisData).color.tickLabel={'10^8','10^9', '10^1^0', '10^1^1', '10^1^2'};
Data(thisData).color.label={'[m^-^3]'};
Data(thisData).y.tick=[70, 80, 90, 100, 120, 140];
Data(thisData).y.tickLabel={'70', '80', '90', '100', '120', '140'};
Data(thisData).y.lim=[60 150];

% Ionosphere Inverse Energy Spectra (Vickrey Model)
thisData=5;
Data(thisData).time = pfisrBeamAvg.energySpectraVickery.time;
Data(thisData).yAxis = pfisrBeamAvg.energySpectraVickery.energyBin;
Data(thisData).zValue = pfisrBeamAvg.energySpectraVickery.energyFlux;
Data(thisData).label = ({'PFISR Inversion','(Vickrey et.al 1982)','Diff. energy flux','[keV]'});
Data(thisData).color.lim=10.^[6 12];
Data(thisData).color.tick=10.^[6 8 10 12];
Data(thisData).color.tickLabel={'10^6','10^8', '10^1^0', '10^1^2'};
Data(thisData).color.label={'[eV/m^2 sr s eV]'};
Data(thisData).y.tick=10.^[3 4 5 5.477];
Data(thisData).y.tickLabel={'1','10','100','300'};
Data(thisData).y.lim=[1*10^3 900*10^3];

% Ionosphere Inverse Energy Spectra (SIC Model)
thisData=6;
Data(thisData).time = pfisrBeamAvg.energySpectraSIC.time';
Data(thisData).yAxis = pfisrBeamAvg.energySpectraSIC.energyBin;
Data(thisData).zValue = pfisrBeamAvg.energySpectraSIC.energyFlux;
Data(thisData).label = ({'PFISR Inversion','SIC Model','Diff. energy flux','[keV]'});
Data(thisData).color.lim=10.^[6 12];
Data(thisData).color.tick=10.^[6 8 10 12];
Data(thisData).color.tickLabel={'10^6','10^8', '10^1^0', '10^1^2'};
Data(thisData).color.label={'[eV/m^2 sr s eV]'};
Data(thisData).y.tick=10.^[3 4 5 5.477];
Data(thisData).y.tickLabel={'1','10','100','300'};
Data(thisData).y.lim=[1*10^3 900*10^3];

% Space Measurements
% B field Z and X from FGM
thisData=7;
Data(thisData).time = thd.magneticVectorFieldFGM.time;
Data(thisData).yAxis = thd.magneticVectorFieldFGM.Bz;
Data(thisData).yAxis2 = thd.magneticVectorFieldFGM.Bx;
Data(thisData).label = ({'Thm-D','B-field','[nT]'});
Data(thisData).y.lim = [min(Data(thisData).yAxis2)-10 max(Data(thisData).yAxis)+5];
Data(thisData).y.tick=[-75, -50, -25, 0, 25, 50, 75];
Data(thisData).y.tickLabel={'-75', '-50', '-25', '0', '25', '50', '75'};
Data(thisData).legend = [{'B_z'},{'B_x'}];

% Energy spectra from thd
thisData=8;
Data(thisData).time = thd.energySpectra.time;
Data(thisData).yAxis = thd.energySpectra.energyBin;
Data(thisData).zValue = thd.energySpectra.energyFlux';
Data(thisData).label = ({'Thm-D','Diff. e^- energy flux','[keV]'});
Data(thisData).color.lim=10.^[6 12];
Data(thisData).color.tick=10.^[6 8 10 12];
Data(thisData).color.tickLabel={'10^6','10^8', '10^1^0', '10^1^2'};
Data(thisData).color.label={'[eV/m^2 sr s eV]'};
Data(thisData).y.tick=10.^[3 4 5 5.477];
Data(thisData).y.tickLabel={'1','10','100','300'};
Data(thisData).y.lim=[1*10^3 900*10^3];

% Ion Anisotropy
thisData=9;
Data(thisData).time = thd.pitchAngleAnisotropy.sst.ionTime;
Data(thisData).yAxis = thd.pitchAngleAnisotropy.sst.ionAnisotropy;
Data(thisData).time2 = thd.pitchAngleAnisotropy.esa.ionTime;
Data(thisData).yAxis2 = thd.pitchAngleAnisotropy.esa.ionAnisotropy;
Data(thisData).label = ({'Thm-D',' i^+ Anisotropy','[a.u.]'});
Data(thisData).y.tick=[-1.0 -0.5 0 0.5 1.0];
Data(thisData).y.tickLabel={'-1.0','-0.5','0','0.5' , '1.0'};
Data(thisData).y.lim = [-1.0 +1.0];
Data(thisData).legend = [{'> 25 keV'}, {'< 25 keV'}];
% Electron Anisotropy
thisData=10;
Data(thisData).time = thd.pitchAngleAnisotropy.sst.electronTime;
Data(thisData).yAxis = thd.pitchAngleAnisotropy.sst.electronAnisotropy;
Data(thisData).time2 = thd.pitchAngleAnisotropy.esa.electronTime;
Data(thisData).yAxis2 = thd.pitchAngleAnisotropy.esa.electronAnisotropy;
Data(thisData).label = ({'Thm-D','e^- Anisotropy','[a.u.]'});
Data(thisData).y.lim = [-1.0 +1.0];
Data(thisData).y.tick=[-1.0, -0.5, 0, 0.5, 1.0];
Data(thisData).y.tickLabel={'-1.0','-0.5','0','0.5' , '1.0'};
Data(thisData).legend = [{'> 30 keV'}, {'< 30 keV'}];

% Wave power 1-1000Hz from FBK
thisData=11;
Data(thisData).time = thd.waveHighFrequencyFB_EFI.time;
Data(thisData).yAxis = thd.waveHighFrequencyFB_EFI.frequency;
Data(thisData).zValue = thd.waveHighFrequencyFB_EFI.power';

% Correcting for 0 values:
minzValue=min(min(((Data(thisData).zValue(Data(thisData).zValue>0)))));
Data(thisData).zValue(Data(thisData).zValue==0)=minzValue;

Data(thisData).label = ({'Thm-D','FB','[Hz]'});
Data(thisData).color.label={'[mV m^-^1 Hz^-^0^.^5]'};
Data(thisData).y.tick=10.^[0 1 2 3];
Data(thisData).y.tickLabel={'1','10','100','10^3'};
Data(thisData).y.lim=[min(Data(thisData).yAxis) max(Data(thisData).yAxis)];
Data(thisData).color.lim=10.^[-2 1];
Data(thisData).color.tick=10.^[-2 -1 0 1];
Data(thisData).color.tickLabel={'10^-^2','10^-^1', '10^0', '10^1'};

% Wave power <1Hz from SCM along X
thisData=12;
Data(thisData).time = thd.waveLowFrequencySCM.time;
Data(thisData).yAxis = thd.waveLowFrequencySCM.frequency;
Data(thisData).zValue = thd.waveLowFrequencySCM.powerGsmX';

% Correcting for 0 values:
minzValue=min(min(((Data(thisData).zValue(Data(thisData).zValue>0)))));
Data(thisData).zValue(Data(thisData).zValue==0)=minzValue;

Data(thisData).label = ({'Thm-D','SCM GSM-X','[Hz]'});
Data(thisData).color.label={'[nT Hz^-^0^.^5]'};
Data(thisData).y.tick=10.^[-2 -1 0];
Data(thisData).y.lim=[min(Data(thisData).yAxis) max(Data(thisData).yAxis)];
Data(thisData).y.tickLabel={'10^-^2','10^-^1','10^0'};
Data(thisData).color.lim=10.^[-4 3];
Data(thisData).color.tick=10.^[-4 -2 0 2];
Data(thisData).color.tickLabel={'10^-^4','10^-^2', '10^0', '10^2'};

% Wave power <1Hz from EFI along Y
thisData=13;
Data(thisData).time = thd.waveLowFrequencyEFI.time;
Data(thisData).yAxis = thd.waveLowFrequencyEFI.frequency;
Data(thisData).zValue = thd.waveLowFrequencyEFI.wavePowerGsmY';

% Correcting for 0 values:
minzValue=min(min(((Data(thisData).zValue(Data(thisData).zValue>0)))));
Data(thisData).zValue(Data(thisData).zValue==0)=minzValue;

Data(thisData).label = ({'Thm-D','EFI GSM-Y','[Hz]'});
Data(thisData).color.label={'[mV m^-^1 Hz^-^0^.^5]'};
Data(thisData).y.tick=10.^[-2 -1 0];
Data(thisData).y.tickLabel={'10^-^2','10^-^1','10^0'};
Data(thisData).y.lim=[min(Data(thisData).yAxis) max(Data(thisData).yAxis)];
Data(thisData).color.lim=10.^[-4 3];
Data(thisData).color.tick=10.^[-4 -2 0 2];
Data(thisData).color.tickLabel={'10^-^4','10^-^2', '10^0', '10^2'};
% Wave power <1Hz from FGM along X
thisData=14;
Data(thisData).time = thd.waveLowFrequencyFGM.time;
Data(thisData).yAxis = thd.waveLowFrequencyFGM.frequency;
Data(thisData).zValue = thd.waveLowFrequencyFGM.powerGsmX';

% Correcting for 0 values:
minzValue=min(min(((Data(thisData).zValue(Data(thisData).zValue>0)))));
Data(thisData).zValue(Data(thisData).zValue==0)=minzValue;

Data(thisData).label = ({'Thm-D','FGM GSM-X','[Hz]'});
Data(thisData).color.label={'[nT Hz^-^0^.^5]'};
Data(thisData).y.tick=10.^[-2 -1 0];
Data(thisData).y.tickLabel={'10^-^2','10^-^1','10^0'};
Data(thisData).y.lim=[min(Data(thisData).yAxis) max(Data(thisData).yAxis)];
Data(thisData).color.lim=10.^[-4 3];
Data(thisData).color.tick=10.^[-4 -2 0 2];
Data(thisData).color.tickLabel={'10^-^4','10^-^2', '10^0', '10^2'};

% Additional Ground Measurements
% GAKO ASI Keogram 
thisData=15;
Data(thisData).time = asi.gako.time;
Data(thisData).yAxis = asi.gako.mlat;
Data(thisData).zValue = asi.gako.intensity';
Data(thisData).label = ({'GAKO 63.05^0N','Keogram','[MLAT]'});
Data(thisData).y.tick=[60, 62, 64, 66, 68, 70];
Data(thisData).y.tickLabel={'60^0N', '62^0N', '64^0N', '66^0N', '68^0N', '70^0N'};
Data(thisData).color.label={'Intensity'};
Data(thisData).y.lim=[60 70];
Data(thisData).color.lim=10.^[3.4 3.8];
Data(thisData).color.tick=10.^[3.4 3.5 3.6 3.7 3.8];
Data(thisData).color.tickLabel={'10^3^.^4','10^3^.^5', '10^3^.^6', '10^3^.^7', '10^3^.^8'};


% Cumulative Energy Fluxes
% PFISR Vickery
thisData=16;
Data(thisData).zValue= diff_to_cumu_flux(...,
    pfisrBeamAvg.energySpectraVickery.energyFlux,...
    pfisrBeamAvg.energySpectraVickery.energyBin);
Data(thisData).zValue=100*Data(thisData).zValue/max(max(Data(thisData).zValue));
Data(thisData).time = pfisrBeamAvg.energySpectraVickery.time;
Data(thisData).yAxis = pfisrBeamAvg.energySpectraVickery.energyBin;

halfEnergyIndex = find_half_energy( Data(thisData).zValue, 1 );
Data(thisData).yValue = Data(thisData).yAxis(halfEnergyIndex)'; % 50% energy flux Energy

Data(thisData).label = ({'PFISR Inversion','(Vickrey et.al 1982)','Norm. cumulative','energy flux','[keV]'});
Data(thisData).color.lim=10.^[0 2];
Data(thisData).color.tick=10.^[0 1 2];
Data(thisData).color.tickLabel={'0','10', '100'};
Data(thisData).color.label={'[%]'};
Data(thisData).y.tick=10.^[3 4 5 5.477];
Data(thisData).y.tickLabel={'1','10','100','300'};
Data(thisData).y.lim=[1*10^3 900*10^3];


% PFISR (SIC Model)
thisData=17;
Data(thisData).time = pfisrBeamAvg.energySpectraSIC.time';
Data(thisData).yAxis = pfisrBeamAvg.energySpectraSIC.energyBin;
Data(thisData).zValue= diff_to_cumu_flux(...,
    pfisrBeamAvg.energySpectraSIC.energyFlux,...
    pfisrBeamAvg.energySpectraSIC.energyBin);
Data(thisData).zValue=100*Data(thisData).zValue/max(max(Data(thisData).zValue));

halfEnergyIndex = find_half_energy( Data(thisData).zValue, 1 );
Data(thisData).yValue = Data(thisData).yAxis(halfEnergyIndex)'; % 50% energy flux Energy

Data(thisData).label = ({'PFISR Inversion','SIC Model','Norm. cumulative','energy flux','[keV]'});
Data(thisData).color.lim=10.^[0 2];
Data(thisData).color.tick=10.^[0 1 2];
Data(thisData).color.tickLabel={'0','10', '100'};
Data(thisData).color.label={'[%]'};
Data(thisData).y.tick=10.^[3 4 5 5.477];
Data(thisData).y.tickLabel={'1','10','100','300'};
Data(thisData).y.lim=[1*10^3 900*10^3];


% THEMIS
thisData=18;
Data(thisData).time = thd.energySpectra.time;
Data(thisData).yAxis = thd.energySpectra.energyBin;
Data(thisData).zValue= diff_to_cumu_flux(...,
    thd.energySpectra.energyFlux',...
    thd.energySpectra.energyBin');
Data(thisData).zValue=100*Data(thisData).zValue/max(max(Data(thisData).zValue));

halfEnergyIndex = find_half_energy( Data(thisData).zValue, 1 );
Data(thisData).yValue = Data(thisData).yAxis(halfEnergyIndex)'; % 50% energy flux Energy

Data(thisData).label = ({'Thm-D','Norm. cumulative','energy flux','[keV]'});
Data(thisData).color.lim=10.^[0 2];
Data(thisData).color.tick=10.^[0 1 2];
Data(thisData).color.tickLabel={'0','10', '100'};
Data(thisData).color.label={'[%]'};
Data(thisData).y.tick=10.^[3 4 5 5.477];
Data(thisData).y.tickLabel={'1','10','100','300'};
Data(thisData).y.lim=[1*10^3 900*10^3];

% Parallel Potential
% by estimating loss cone angle flux @ ~1.5deg in plasma sheet 
thisData=19;
Data(thisData).time = parallel_potential.time;
Data(thisData).yAxis = parallel_potential.V;
temp=interp1(Data(16).time, Data(16).yValue, Data(18).time);
Data(thisData).yAxis2 = (temp-Data(18).yValue)/1000;
Data(thisData).label = ({'\Delta \phi_|_|','[kV]'});
Data(thisData).y.tick=[-20, 0, 20, 40, 80];
Data(thisData).y.tickLabel={'-20','0', '20','40','80'};
Data(thisData).y.lim = [-20 +80];
Data(thisData).legend = [{'Difference of Peaks'},{'Difference of Median Energy'}];


end

