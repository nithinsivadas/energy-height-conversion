function [ rootPathStr, dataPathStr, computer ] = initialize_root_path()
%% initialize_root_path.m Generates the root git-hub path based on the computer
%--------------------------------------------------------------------------
% Output
%-------
% rootPathStr - A string containg git directory path in the present computer
%----------------------------------------------------------------------------
% Modified: 
% Created : 10th Feb 2017
% Author  : Nithin Sivadas
% Ref     :
% Update  : v1.1 10/31/2017 Added Nithin-Carbon, Also changed the
% comparison function to strcmp, instead of just "=="
%----------------------------------------------------------------------------
computer=get_computer_name();
if nargin<1
    if strcmp(computer,'nithin-surface')
        rootPathStr='C:\Users\Nithin\Documents\GitHub\';
        dataPathStr = 'C:\Users\nithin\Documents\GitHub\LargeFiles';
    elseif strcmp(computer,'nithin-carbon')
        rootPathStr='C:\Users\nithin\Documents\GitHub\';
        dataPathStr = 'C:\Users\nithin\Documents\GitHub\LargeFiles';
    elseif strcmp(computer,'aurora1-optiplex-780')
        rootPathStr='/home/nithin/Documents/git-repos/';
        dataPathStr = '/home/nithin/Documents/git-repos/Largefiles/';
    elseif strcmp(computer,'scc-lite')
        rootPathStr = '/usr3/graduate/nithin/git-repos/';
        dataPathStr = '/scratch/nithin/Data/';
    else
        error('Computer not initialized to run energy-height-conversion/Tools');
    end
end

end

