%% Generate omni.h5 (Linux)

h5FileStr = '/home/nithin/Documents/git-repos/LargeFiles/omni/omni.h5';
localStorePath = [initialize_root_path,'LargeFiles',filesep,'omni',filesep,'ASC'];
setCalculateGW = true;
tic
[status] = create_omni_HDF5(h5FileStr, localStorePath, setCalculateGW);
toc

%% Check omni

[maginput,time,header,units] = generate_maginput(h5FileStr,'26-Mar-2008 11:00','26-Mar-2008 12:00');