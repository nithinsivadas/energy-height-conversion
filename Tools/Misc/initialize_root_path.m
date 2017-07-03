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
%----------------------------------------------------------------------------
computer=getenv('computername');
if nargin<1
    if computer=='NITHIN-SURFACE'
        rootPathStr=['C:\Users\Nithin\Documents\GitHub\'];    
    else
        rootPathStr=['/home/nithin/Documents/git-repos/'];
    end
end;

end

