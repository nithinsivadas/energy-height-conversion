%% POES, NOAA-15, satellite trajectory for a year

folderStr = 'G:\My Drive\Nithin\Job Search\JEPostDoc\Data\POES\';
fileStr = get_files_in_folder(folderStr,'*.bin');

%%
orbit.lat = [];
orbit.lon = [];
orbit.alt = [];
orbit.time = [];
orbit.mep0E1 = [];
orbit.mep0E2 = [];
orbit.mep0E3 = [];

multiWaitbar('Compiling POES data...',0);
n = length(fileStr);
dn = 1./n;
for i = 1:1:n
    iFile = [folderStr,fileStr{i}];
    [poes, Table, outputFile] = POES_extract_binary_faster(iFile);
    orbit.time = [orbit.time; poes.time];
    orbit.lat = [orbit.lat; poes.lat];
    orbit.lon = [orbit.lon; poes.lon];
    orbit.alt = [orbit.alt; poes.alt];
    orbit.mep0E1 = [orbit.mep0E1; poes.mep0E1];
    orbit.mep0E2 = [orbit.mep0E2; poes.mep0E2];
    orbit.mep0E3 = [orbit.mep0E3; poes.mep0E3];
    multiWaitbar('Compiling POES data...','Increment',dn);
end
multiWaitbar('Close all');