function dlFITS(myurl,final_dir,times,varargin)
% dlFITS
% by John Swoboda
% This function will find the data within a certain time period and down it
% to specific directory.
%% Inputs
% myurl - The url that all the files are located.
% final_dir - The directory that the data will be downloaded to.
% times - A cell array of date strings that will hold the time limits
% wl - The desired wavelength of the optical data.
%% Example
% myurl = 'http://amisr.asf.alaska.edu/PKR/DASC/RAW/2012/20121124/';
% final_dir =  '/Volumes/Research/eng_research_irs/PINOT/Data_Fusion/FusedData/';
% times = {'11/24/2012 6:00:00','11/24/2012 6:15:00'};
% wl = 558;
% dlFITS(myurl,final_dir,times,wl)
%%
p = inputParser;
addOptional(p,'wl',[])
p.parse(varargin{:})
U = p.Results;

if ~exist(final_dir,'dir')
    error(['Your output directory ',final_dir,' does not exist'])
end
disp(['outputting files to ',final_dir])
%%
allfiles = htmlfindfile(myurl,'*\.FITS');

red_file_list = fitslistparce(allfiles,times,U.wl);

nfile = length(red_file_list);
if nfile > 100
    warning(['Attempting to download ',int2str(nfile),' files, this may take a long time and use a lot of Hard drive space.'])
end

for k = 1:nfile
    temp_filename = [myurl,red_file_list{k}];
    temp_fileput = fullfile(final_dir,red_file_list{k});

    updatestr = [red_file_list{k},' ', int2str(k),' / ',int2str(nfile)];
    if exist(temp_fileput,'file')
        disp(['skipping already existing ',updatestr])
        continue
    else
        disp(['downloading ',updatestr])
    end

    try
        websave(temp_fileput,temp_filename);
    catch
        urlwrite(temp_filename,temp_fileput);
    end
end %for

end %function
