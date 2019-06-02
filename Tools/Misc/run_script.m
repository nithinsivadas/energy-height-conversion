function run_script(scriptFnHandle)
%% run_script.m Runs a script written in the form of a function
% This function can be run from the terminal. It records relevant time
% , errors encountered while executing the scriptFnHandle that is passed,
% and sends an email to a preset email id after the function is executed.  
% INPUT:
%       scriptFnHandle - Handle to a script function
%                        e.g. @AMISR_exp_script 
%--------------------------------------------------------------------------
% Created by: Nithin Sivadas
% Created on: 2nd June 2019
% Updated on: NA
%--------------------------------------------------------------------------

if exist('record.txt','file')
    if ispc
        system('del record.txt')
    else
        system('rm record.txt');
    end
end

diary record.txt
errorFlag = 0;
startTimeStr = (datevec(datetime('now')));
disp(['Start time: ',datestr(startTimeStr), 10]);
if exist('scriptFnHandle','var')==1
    funcStr=func2str(scriptFnHandle);
    try 
        %% Script 
        scriptFnHandle();

    catch ME
        disp('-----------------------------------------');
        disp('            Encountered Error');
        disp('-----------------------------------------');
        errorFlag = 1;
       getReport(ME)
    end
else
    disp([10 'Error: No input function handle' 10]);
    errorFlag = 1;
    funcStr = 'No function handle, script has';
end
endTimeStr = (datevec(datetime('now')));
disp([10,'End time: ',datestr(endTimeStr),10]);
disp('-----------------------------------------');
disp('              Complete');
disp('-----------------------------------------');
diary off

send_email('attachments','record.txt','scriptName',...
   funcStr,'errorFlag',errorFlag);

if ispc
    system('del record.txt')
else
    system('rm record.txt');
end

end