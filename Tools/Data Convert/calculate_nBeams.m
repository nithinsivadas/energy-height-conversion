function [ nBeams ] = calculate_nBeams(pfisrGD, coordinateNo)
%calculate_nBeams Calculate the number of beams in PFISR experiment
% using the range, or slant_height coordinate - which is same for all beams
%   Detailed explanation goes here
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

