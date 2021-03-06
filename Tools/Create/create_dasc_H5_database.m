function [status, err] = create_dasc_H5_database...
    (timeMinStr,timeMaxStr,outputH5File,writeModeStr)
%create_amisr_web_experiment_H5_database.h5 Downloads from DASC FTP the data
% regarding DASC experiments, and stores it in an HDF5 file
%-------------------------------------------------------------------------
% Input
%------
%       timeMinStr   - Start time string
%       timeMaxStr   - Stop time string
%       outputH5File - Output database file name
%       writeModeStr -  NA
%----------------------------------------------------------------------------
% Output
%---------
%       status - 'Success' or 'Failed'
%----------------------------------------------------------------------------
% Modified: 15th Aug 2019, 
%           2nd Dec 2020 - changed hdf5 variable wavelength [dependent
%           functions will need to be also altered]
% Created : Unknown
% Author  : Nithin Sivadas
% Ref     :
% Comments:
%
%----------------------------------------------------------------------------

if nargin<3
    outputH5File = 'dascDatabase.h5';
end

if nargin<2 || isempty(timeMaxStr)
    timeMaxStr = '20 Jan 2018';
end
if nargin <1 || isempty(timeMinStr)
    timeMinStr = '31 Dec 2017';
end

status = 'failed';
err = struct();

if isfile(outputH5File)
   movefile(outputH5File, [outputH5File,'_bak']);
end

date1 = datetime(datenum(timeMinStr),'ConvertFrom','datenum');
date2 = datetime(datenum(timeMaxStr),'ConvertFrom','datenum');

% Can be edited to include more cameras
instrumentStr = 'DASC';

yearArr = year(date1):year(date2);
for iYear = 1:1:length(yearArr)

   try 
        disp(['Year: ',num2str(yearArr(iYear))]); %Marker
        
        if iYear == 1
            timeMinStr1 = timeMinStr;
        else
            timeMinStr1 = ['01 Jan ',num2str(yearArr(iYear))];
        end
        
        if iYear == length(yearArr)
            timeMaxStr1 = timeMaxStr;
        else
            timeMaxStr1 = ['31 Dec ',num2str(yearArr(iYear)),' 23:59:59.999'];
        end

        [data, err] = get_DASC_FITS_times(timeMinStr1, timeMaxStr1);
        
        if ~isempty(data.time)
        
            disp('Writing to H5'); %Marker
            data.time1=posixtime(data.time);
            
            waveUnique = unique(data.wavelength);
            waveID = zeros(size(data.wavelength));
            waveCode = string();
            for i=1:size(waveUnique,1)
            waveID = waveID +...
                strcmp(data.wavelength,waveUnique(i)).*...
                find(strcmp(waveUnique,waveUnique(i)));
            waveCode(i,1)=waveUnique(i);
            waveCode(i,2)=num2str(i);
            end
            
            h5createwritestr(outputH5File, ['/',instrumentStr,'/',num2str(yearArr(iYear)),...
            '/wavelengthCode'],  (waveCode)');
            h5writeatt(outputH5File, ['/',instrumentStr,'/',num2str(yearArr(iYear)),...
            '/wavelengthCode'], 'Description', ['2D array of unique image wavelengths,',...
            'recorded during this day, with its corresponding integer identifier.']);

            write_h5_dataset(outputH5File,['/',instrumentStr,'/',num2str(yearArr(iYear)),...
            '/wavelength'],waveID,1,true);
            h5writeatt(outputH5File, ['/',instrumentStr,'/',num2str(yearArr(iYear)),...
            '/wavelength'], 'Description', ['Wavelength of each recorded image, ',...
            'identified using an integer code defined in /wavelengthCode.']);
        
            write_h5_dataset(outputH5File,['/',instrumentStr,'/',num2str(yearArr(iYear)),...
            '/time'],(data.time1),1,true);
            h5writeatt(outputH5File,['/',instrumentStr,'/',num2str(yearArr(iYear)),...
            '/time'],'Description','[1 x nT] Time in POSIX units');       
        end
        
   catch ME
       
        disp([[10 'Error in documenting Year: '],num2str(yearArr(iYear))]);
        disp(getReport(ME));
        
   end
end
        
        
status = 'Success';

end
