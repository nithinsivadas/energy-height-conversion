function create_dasch5_script()
    % Generate DASC database from ftp://optics.gi.alaska.edu/
    % Store the data in root data path
    % run this from linux commad line as follows:
    % runMatlabScript @create_dasch5_script
    [~, dataPath] = initialize_root_path();
%     dataPathStr = "/media/nithin/PFISR_002_006/Nithin/";
    outputFileStr = 'dascDatabase.h5';
%     timeMinStr = '01 Jan 2007';
%     timeMaxStr = '31 Dec 2019';
    timeMinStr = '01 Jan 2007';
    timeMaxStr = '31 Dec 2019';
    [status, error] = create_dasc_H5_database(timeMinStr, timeMaxStr, [dataPath,outputFileStr]);
    disp([10 'Status: ']);
    cellstr(status)
    disp([10 'Data stored in ',dataPath,outputFileStr, 10]);
end