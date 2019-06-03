function AMISR_exp_script()
    % Generate AMISR experiment database from amisr.com
    % Store the data in root data path
    [~, dataPath] = initialize_root_path();
%     dataPathStr = "/media/nithin/PFISR_002_006/Nithin/";
    outputFileStr = 'amisrWebDatabase.h5';
    timeMinStr = '1 Dec 2006';
    timeMaxStr = '1 Jun 2019';
    create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,[dataPath,outputFileStr]);
    disp([10 'Data stored in ',dataPath,outputFileStr, 10]);
end