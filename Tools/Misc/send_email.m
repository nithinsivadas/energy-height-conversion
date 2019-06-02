function send_email(varargin)
%send_email Sends an email from your gmail id set by preferences. 
%-------
% INPUT:
%       to          - Recepient's email id
%       subject     - Subject of the email
%       message     - Message to be sent [has a default format]
%       scriptName  - Name of the script executing this function, will be
%                     included in the default message
%       attachments - File paths to be attached, can be a cell array of
%                     strings

% Make sure your preferences are already set
%setpref('Internet','SMTP_Server','smtp.gmail.com');
%setpref('Internet','E_mail','nithin.iitm@gmail.com');
%setpref('Internet','SMTP_Username','myusername');
%setpref('Internet','SMTP_Password','mypassword');

p = inputParser;

addParameter(p,'to','nithin@bu.edu');
addParameter(p,'subject','Code executed');
addParameter(p,'attachments',nan);
addParameter(p,'errorFlag',0);
addParameter(p,'scriptName','A MATLAB script');
nowTimeStr = datestr(clock);

addParameter(p,'message',nan);

parse(p,varargin{:});

scriptName = p.Results.scriptName;

if isnan(p.Results.message)
    if p.Results.errorFlag == 0
    message = ['Hello'...
        10 10 scriptName,' has been executed'...
        10 'on ',get_computer_name...
        10 'at ',nowTimeStr...
        10 10 'Thank you'...
        10 '- MATLAB ',version];
        subject = ['[MATLAB] ',p.Results.subject,' ',nowTimeStr];
    else
        message = ['Hello'...
        10 10 scriptName,' has failed to execute'...
        10 'on ',get_computer_name...
        10 'at ',nowTimeStr...
        10 10 'Sorry!'...
        10 '- MATLAB ',version];
        subject = ['[MATLAB] ','Error encountered',' ',nowTimeStr];
    end
end

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

if ~isnan(p.Results.attachments)
    sendmail(p.Results.to,subject,message,p.Results.attachments);
else
    sendmail(p.Results.to,subject,message);
end


end