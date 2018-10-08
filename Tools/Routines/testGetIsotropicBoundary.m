%% testGetIsotropicBoundary.m
clear all;

%% Initialize
magFieldModelNo = 7;
homeDir = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\';
inputH5FileStr=[homeDir,'20080326.001_bc_15sec-energyFlux.h5'];
timeStr = '26 Mar 2008 11:10';
energy = 100; % [keV]
%%

magFieldModelStr=find_irbem_magFieldModelStr(magFieldModelNo);
magData = get_2D_plot_inputs_time_independent(inputH5FileStr,...
    'plotModeStr','MagneticFieldMap','magFieldModelStr',magFieldModelStr);
thisTimeIndx=find_time(magData.time,timeStr);
magData.BgeoEq = readh5_variable_at_time(inputH5FileStr,...
            'BgeoEq',['/magneticMap/',magFieldModelStr,'/'],thisTimeIndx)';
magData.BmagEq = readh5_variable_at_time(inputH5FileStr,...
            'BmagEq',['/magneticMap/',magFieldModelStr,'/'],thisTimeIndx);
magData.gradBmagEq = readh5_variable_at_time(inputH5FileStr,...
            'gradBmagEq',['/magneticMap/',magFieldModelStr,'/'],thisTimeIndx)';
magData.diffBEq = permute(readh5_variable_at_time(inputH5FileStr,...
            'diffBEq',['/magneticMap/',magFieldModelStr,'/'],thisTimeIndx),[3,2,1]);
        
%%
[Kc] = get_isotropic_boundary(magData.BgeoEq,magData.BmagEq,...
    magData.gradBmagEq,magData.diffBEq,energy);