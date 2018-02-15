function [fileStr] = get_files_in_folder(dirName,wildcard)
%get_files_in_folder.m Produces a cell array of all file names in folder 
%   Detailed explanation goes here
if nargin<2
    wildcard='*';
end

files = dir([dirName,'/',wildcard]);
fileIndex = find(~[files.isdir]);
fileTempCells  = struct2cell(files);
fileStr = fileTempCells(1,fileIndex);

end

