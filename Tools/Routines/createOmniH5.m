%% Generate omni.h5 (Linux)

h5FileStr = '/home/nithin/Documents/git-repos/LargeFiles/omni.h5';
localStorePath = [initialize_root_path,'LargeFiles',filesep,'omni',filesep,'ASC'];
setCalculateGW = true;
[data,status] = create_omni_HDF5(h5FileStr, localStorePath, setCalculateGW);