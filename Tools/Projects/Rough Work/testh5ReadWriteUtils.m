%% Routine building up some reliable h5read and h5write tools
% fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20100528.001_bc_2min-energyFlux.h5';
% fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\20080326.001_bc_15sec-energyFlux.h5';
% %% Test existing h5read_data.m
% % [data, h5DataLoc] =  h5read_data(fileStr);
% 
% %% Writing function to find the address of all Datasets
% 
% % Returns a string array of all dataset paths in HDF5 file
% datasetPaths = read_h5_dataset_paths(fileStr,'/');
% 
% % Reads all data within a specified group (or all data) in HDF5 file
% data=read_h5_data(fileStr,'/magneticMap/TS96');
% 
% % Reading h5 variable at particular time

%%
fileStr = 'G:\My Drive\Research\Projects\Paper 2\Data\EEA List\test.h5';
% time = (datenum('26-Mar-2008 13:11'):(1/3600):datenum('26-Mar-2008 13:20'))';
% varValue = 2*ones(length(time),26,192);
% 
% write_h5_dataset(fileStr, '/testData3/time', time,1,'true','true');
% write_h5_dataset(fileStr, '/testData3/Ne', varValue,1,'true','true');
write_h5_dataset_attribute(fileStr, '/testData3/time', [],[],[], 'time');
write_h5_dataset_attribute(fileStr, '/testData3/Ne',...
    'Test electron denisty','nTime x nBeam x nAlt','m-3');