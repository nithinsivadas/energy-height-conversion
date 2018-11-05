%% Add thg to HDF5

%% Initialization
% outputFileStr = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20101018.001_bc_15sec-energyFlux_v85km.h5';
outputFileStr ='/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-energyFlux_v85.h5';
siteName = 'gako';

%% Execute
status = create_thg_hdf5(siteName,outputFileStr);
multiWaitbar('Close All');