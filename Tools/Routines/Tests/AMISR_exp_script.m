function AMISR_exp_script()
    dataStoreDir = "/media/nithin/PFISR_002_006/Nithin/";
    outputFileStr = "amisrWebDatabase.h5";
    timeMinStr = '1 Dec 2006';
    timeMaxStr = '1 Jun 2019';
    create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,strcat(dataStoreDir,outputFileStr));
end