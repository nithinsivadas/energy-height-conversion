%% Linux create omni script
diary record.txt
startTimeStr = (datevec(datetime('now')));
disp(['Start time: ',datestr(startTimeStr)]);
try 
    %% Your Code
    dataStoreDir = "/media/nithin/PFISR_002_006/Nithin/";
    outputFileStr = "amisrWebDatabase.h5";
    timeMinStr = '1 Dec 2006';
    timeMaxStr = '1 Jun 2019';
    status = create_amisr_web_experiment_H5_database(timeMinStr,timeMaxStr,61,strcat(dataStoreDir,outputFileStr));
    % Code ends
catch ME
    disp('-----------------------------------------');
    disp('            Encountered Error');
    disp('-----------------------------------------');
    ME
end
endTimeStr = (datevec(datetime('now')));
disp(['End time: ',datestr(endTimeStr)]);
disp('-----------------------------------------');
disp('              Complete');
disp('-----------------------------------------');
diary off
send_email('attachments','record.txt','scriptName',mfilename);
system('rm record.txt');
exit