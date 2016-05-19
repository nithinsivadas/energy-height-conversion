function [data] = combine_energies(data_esa,data_sst)
%combine_energies Combines the esa and sst electron energy spectras
%  Inputs:
%  data_esa and data_sst which both contain 
%               E[NxM]: energy flux [eV /cm^2 sec sr  eV]
%           time [Nx1]: time in matlab units
%            ebin[1xM]: energy bins [1xM]
%  Outputs:
%  data
%  data.E[N x (M1+M2)]: the combined matrix of esa and sst energies
%  data.time[Nx1]     : the time instants in sst data set
%  data.ebin[1x M1+M2]: the combined energy bin array
% 
% Note that we remove one inconsistent energy bin from esa (~30 keV) so as
% not to overlap with the sst low energy value (~30 keV)

% 12 May 2016
% Nithin Sivadas


ebin=[fliplr(data_esa.ebin(2:end)),data_sst.ebin];
E   =[fliplr(data_esa.E(:,2:end)),data_sst.E];
if isequal(data_esa.time,data_sst.time)~=1
    display('Time axis is not exactly equal for sst and esa');
end
data.E=E;
data.ebin=ebin;
data.time=data_sst.time;
end

