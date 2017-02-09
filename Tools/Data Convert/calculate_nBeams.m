function [ nBeams ] = calculate_nBeams(pfisrGD, coordinateNo)
%% calculate_nBeams Calculate the number of beams in PFISR experiment
% using the range, or slant_height coordinate - which is same for all beams
%--------------------------------------------------------------------------
% Input
%------
% pfisrGD      : Geodata object of PFISR
% coordinateNo : The coordinate number to be used to estimate the numbero
%                of beams (usually the slant height)
%--------------------------------------------------------------------------
% Output
%-------
% nBeams : The number of beams present in the experiment by PFISR
%----------------------------------------------------------------------------
% Modified: 24th Jan 2017 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------

    if nargin<2
        coordinateNo = 1;
    end;

    temp = pfisrGD.dataloc(1,coordinateNo);
    j=1;
    while (pfisrGD.dataloc(j+1,coordinateNo)==temp)
        j=j+1;
    end;    
    nBeams=j; % Total number of beams

end

