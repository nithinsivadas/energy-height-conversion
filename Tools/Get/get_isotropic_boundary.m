function [Kc] = get_isotropic_boundary(Bgeo,Bmag,...
    gradBmag,diffB,energy)
%GET_ISOTROPIC_BOUNDARY Gets the ratio of the Rcurve to Gyroraidus of electron, 
% which if ~ 8 defines the isotropic boundary

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

