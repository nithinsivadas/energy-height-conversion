function create_dasch5_script()
    % Generate DASC database from ftp://optics.gi.alaska.edu/
    % Store the data in root data path
    % run this from linux commad line as follows:
    % runMatlabScript @create_dasch5_script
    [~, dataPath] = initialize_root_path();
%    dataPath = "/projectnb/semetergrp/nithin/Data/";
    outputFileStr = 'dascDatabase.h5';

%     timeMinStr = '1 Jan 2007';
%     timeMaxStr = '1 Jan 2008';

    timeMinStr = '01 Jan 2007';
    timeMaxStr = '29 Nov 2020 11:59:59.999';
    
    [status, err] = create_dasc_H5_database(timeMinStr, timeMaxStr, [dataPath,outputFileStr]);
    disp([10 'Status: ']);
    cellstr(status)
    
    disp([10 'Errors: ']);
    disp([err.directory.name,err.directory.message]);
    disp([err.subdirectory.name,err.subdirectory.message]);
    
    disp([10 'Data stored in ',dataPath,outputFileStr, 10]);
end