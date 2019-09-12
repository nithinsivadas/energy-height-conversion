%% Download DASC days for calibration

storeDir = 'G:\My Drive\Research\Projects\Paper 3\Data\WorkDir\';

timeStr{1} = '14 Dec 2010';
timeStr{2} = '12 Nov 2016';

download_DASC_FITS(timeStr{1},storeDir);
fileStr1=find_quite_time_frame(timeStr{1}, storeDir);
fileStr1

download_DASC_FITS(timeStr{2},storeDir);
fileStr2=find_quite_time_frame(timeStr{1}, storeDir);
fileStr2




function fileStr = find_quite_time_frame(timeStr, storeDir)

    dirName = strcat(storeDir,datestr(timeStr,'yyyymmdd'));
    [fileStrList] = get_files_in_folder(dirName,'*.FITS');
    fullFileStr = string(strcat(dirName,filesep,fileStrList)');
    
    multiWaitbar('Calculating ...', 0');
    di = 1./length(fullFileStr);
    
    for i = 1:1:length(fullFileStr)
        multiWaitbar('Calculating ...','Increment',di);
        tempData = fitsread(fullFileStr(i));
        intensityArr(i) = nansum(tempData(:));
    end
    multiWaitbar('Calculating ...',1);
    
    [~,I] = min(intensityArr);
    fileStr = fullFileStr(I);
    
end