%% Test _ email


% startTimeStr = (datevec(datetime('now')));
% disp(['Start time: ',datestr(startTimeStr),10]);

disp('We are testing send_email function');

disp('Checking the internet preferences');

disp(['SMTP_Server   :',getpref('Internet','SMTP_Server')]);
disp(['Email         :',getpref('Internet','E_mail')]);
disp(['SMTP_Username :',getpref('Internet','SMTP_Username')]);
disp('If this much is correct, then the password should be stored fine.');

disp('Displaying the current directory');
disp(pwd);

endTimeStr = (datevec(datetime('now')));
% disp(['Start time: ',datestr(endTimeStr),10]);


