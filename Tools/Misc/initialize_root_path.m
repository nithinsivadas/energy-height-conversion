function [ rootPathStr ] = initialize_root_path()
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
computer=getenv('computername');
if nargin<1
    if strcmp(computer,'NITHIN-SURFACE')
        rootPathStr=['C:\Users\Nithin\Documents\GitHub\'];    
    elseif strcmp(computer,'NITHIN-CARBON')
        rootPathStr=['C:\Users\nithin\Documents\GitHub\'];    
    else 
        rootPathStr=['/home/nithin/Documents/git-repos/'];
    end
end

end
