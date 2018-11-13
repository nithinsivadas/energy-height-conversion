%% Add thg to HDF5

%% Initialization
% outputFileStr = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20101018.001_bc_15sec-energyFlux_v85km.h5';
sites={'fykn','inuv','whit','mcgr','kian'};
multiWaitbar('Processing different cameras',0);
for i=1:1:length(sites)
multiWaitbar('Processing different cameras','Increment',1./length(sites));
outputFileStr ='/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-energyFlux_v85.h5';
siteName = sites{i};
disp(['Processing ',upper(siteName)]);
status = create_thg_hdf5(siteName,outputFileStr);
end
multiWaitbar('Close All');