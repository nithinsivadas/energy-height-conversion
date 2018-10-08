function [Kc] = get_isotropic_boundary(Bgeo,Bmag,...
    gradBmag,diffB,energy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Calculating G
C = define_universal_constants; 
m = C.me;
E = energy*C.e*10^3; %J
V = (2.*E./m).^0.5;
G = (m.*V./(C.e*Bgeo(3)*10^-9))./C.RE;

[~,~,~,~,Rcurv,~,~,~] =...
    onera_desp_lib_compute_grad_curv_curl(Bgeo,Bmag,gradBmag,diffB);

Kc = (Rcurv)./G;
end

