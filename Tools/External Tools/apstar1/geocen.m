function [pos2, tlight] = geocen (pos1, pe)

% this function moves the origin of coordinates from the
% barycenter of the solar system to the center of mass of the
% earth, i.e., this function corrects for parallax.

% input

%  pos1 = position vector, referred to origin at solar system
%         barycenter (au)
%  pe   = position vector of center of mass of the earth,
%         referred to origin at solar system barycenter (au)

% output

%   pos2   = position vector, referred to origin at center of
%            mass of the earth (au)
%   tlight = light time from body to earth in days

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 173.1446333;

for j = 1:1:3
   pos2(j) = pos1(j) - pe(j);
end

tlight = sqrt(pos2(1)^2 + pos2(2)^2 + pos2(3)^2) / c;

