function create_omni_script()
    % Generate AMISR experiment database from amisr.com
    % Store the data in root data path
    % run this from linux commad line as follows:
    % runMatlabScript @AMISR_exp_script
    [~, dataPath] = initialize_root_path();
%     dataPathStr = "/media/nithin/PFISR_002_006/Nithin/";
    outputFileStr = 'omni.h5';
        
    status = create_omni_HDF5([dataPath,outputFileStr],[dataPath,'omni',filesep,'ASC']);
    disp([10 'Status: ']);
    cellstatus(status);
    disp([10 'Data stored in ',dataPath,outputFileStr, 10]);
end