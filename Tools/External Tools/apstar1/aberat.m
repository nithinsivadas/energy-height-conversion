function pos2 = aberat (pos1, ve, tlight)

% this function corrects position vector for aberration of light.
% algorithm includes relativistic terms.  see murray (1981)
% mon. notices royal ast. society 195, 639-648.

% input

%  pos1   = position vector, referred to origin at center of
%           mass of the earth, components in au
%  ve     = velocity vector of center of mass of the earth,
%           referred to origin at solar system barycenter,
%           components in au/day
%  tlight = light time from body to earth in days
%           if tlight = 0, this function will compute tlight

% output

%   pos2 = position vector, referred to origin at center of
%          mass of the earth, corrected for aberration,
%          components in au

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos2 = zeros(3, 1);

c = 173.1446333;

tl = tlight;

p1mag = tl * c;

if (tl ~= 0)
   vemag = sqrt(ve(1)^2 + ve(2)^2 + ve(3)^2);
else
   p1mag = sqrt(pos1(1)^2 + pos1(2)^2 + pos1(3)^2);

   tl = p1mag / c;

   vemag = sqrt(ve(1)^2 + ve(2)^2 + ve(3)^2);
end

beta = vemag / c;

pdotv = pos1(1) * ve(1) + pos1(2) * ve(2) + pos1(3) * ve(3);

cosd = pdotv / (p1mag * vemag);

gammai = sqrt(1 - beta^2);

p = beta * cosd;

q = (1 + p / (1 + gammai)) * tl;

r = 1 + p;

for j = 1:1:3
   pos2(j) = (gammai * pos1(j) + q * ve(j)) / r;
end

