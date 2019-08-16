function [status, error] = create_dasc_H5_database...
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
% Modified: 15th Aug 2019
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

[data, error] = get_DASC_FITS_times(timeMinStr, timeMaxStr);

% Can be edited to include more cameras
instrumentStr = 'DASC';

data.time=posixtime(data.time);

try h5create(outputH5File,['/',instrumentStr,'/time'],size(data.time'),'ChunkSize',[1 80],'Deflate',9);catch ME; end
try h5create(outputH5File,['/',instrumentStr,'/wavelength'],size(data.wavelength'),'ChunkSize',[1 80],'Deflate',9);catch ME; end

h5write(outputH5File,['/',instrumentStr,'/time'],(data.time'));
h5write(outputH5File,['/',instrumentStr,'/wavelength'],(data.wavelength)');
h5_create_writestr(outputH5File,['/',instrumentStr,'/url'],cellstr(data.url)');


% try h5create(outputH5File,['/',instrumentStr,'/file'],size(data.file'),'ChunkSize',[1 80],'Deflate',9);catch ME; end
% try h5create(outputH5File,['/',instrumentStr,'/date'],size(data.date'),'ChunkSize',[1 80],'Deflate',9);catch ME; end
% try h5create(outputH5File,['/',instrumentStr,'/url'],size(data.url'),'ChunkSize',[1 80],'Deflate',9);catch ME; end

% hdf5write(outputH5File,['/',instrumentStr,'/url'],(data.url)','WriteMode','append','Deflate',9);

% ['/',instrumentStr,'/file'],(data.file)',...
% ['/',instrumentStr,'/date'],(data.date)',...


status = 'Success';

end
