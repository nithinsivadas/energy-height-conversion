function pos2 = sunfld (pos1, pe)

% this function corrects position vector for the deflection
% of light in the gravitational field of the sun. see murray (1981)
% mon. notices royal ast. society 195, 639-648. this function valid
% for bodies within the solar system as well as for stars.

% input

%  pos1 = position vector, referred to origin at center of mass
%         of the earth, components in au
%  pe   = position vector of center of mass of the earth,
%         referred to origin at center of mass of the sun,
%         components in au

% output

%  pos2 = position vector, referred to origin at center of mass
%         of the earth, corrected for gravitational deflection,
%         components in au

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mau = 1.49597870e11;

gs = 1.32712438e20;

c = 299792458.0;

% construct vector pq between sun and body

for j=1:1:3
    pq(j) = pe(j) + pos1(j);
end

% compute vector magnitudes and unit vectors

pmag = sqrt (pos1(1)^2 + pos1(2)^2 + pos1(3)^2);
emag = sqrt (pe(1)^2 + pe(2)^2 + pe(3)^2);
qmag = sqrt (pq(1)^2 + pq(2)^2 + pq(3)^2);

for j=1:1:3
    phat(j) = pos1(j) / pmag;
    ehat(j) = pe(j) / emag;
    qhat(j) = pq(j) / qmag;
end

% compute dot products of vectors

pdotq = phat(1) * qhat(1) + phat(2) * qhat(2) + phat(3) * qhat(3);
edotp = ehat(1) * phat(1) + ehat(2) * phat(2) + ehat(3) * phat(3);
qdote = qhat(1) * ehat(1) + qhat(2) * ehat(2) + qhat(3) * ehat(3);

% compute scalar factors

fac1 = 2 * gs / (c * c * emag * mau);

fac2 = 1 + qdote;

% construct corrected position vector pos2

for j=1:1:3
    p2j = phat(j) + fac1 * (pdotq * ehat(j) - edotp * qhat(j)) / fac2;
    
    pos2(j) = p2j * pmag;
end

