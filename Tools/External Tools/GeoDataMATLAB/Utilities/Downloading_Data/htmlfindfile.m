function file_list = htmlfindfile(url,file2find)
% htmlfindfile.m
% file_list = htmlfindfile(url,file2find)
% by John Swoboda
% This function will parce through the html string that is in the url
% address.  The files names will be returned as strings in the cell
% file_list that the user can more easily parce and then use with urlwrite
% to get the files downloaded.  The file2find string is being put into a 
% regular expression so one needs to adjust the string accordingly.  
%% Inputs
% url - a string that contains the URL of the website that one wants to
% look at.
% file2find - A string that will be used as a part in the regular
% expression to find the desired file(s).  Use '*\.FITS' for example to find
% FITS files.
%% Outputs 
% file_list - A Nx1 cell array of strings with the file names.
%% example 
% url = 'http://amisr.asf.alaska.edu/PKR/DASC/RAW/2012/20121124/'
% FITS_file_list = htmlfindfile(url,'*\.FITS');
%% Reference 
% http://stackoverflow.com/questions/11126721/using-matlab-to-parse-html-for-url-in-anchors-help-fast

my_tok = ['<a href="([^"]',file2find,')">'];
html = urlread(url);

file_list = regexp(html,my_tok, 'tokens');
file_list = [file_list{:}]';

end