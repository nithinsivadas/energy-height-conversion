function [data] = read_tplot_ascii(filename1, filename2 )
%read_tplot Reads the ascii files and outputs data structure which contain 
%           enregy, time, and flux
%   Inputs:
%          filename1 - contains the energy flux files
%          filename2 - contains the vertical energy axis
%    Outputs:
%               data.E[NxM]: energy flux [eV /cm^2 sec sr  eV]
%           data.time [Nx1]: time in matlab units
%            data.ebin[1xM]: energy bins [1xM]
% 12 May 2016
% Nithin Sivadas

A = importdata(filename1,' ');
V = importdata(filename2,' ');

E = A.data(:,:);
time = datenum(A.rowheaders);
ebin = V(1,:);

data.time = time;
data.ebin = ebin;
data.E    = E;

end

