% example 4: Besselian and Julian epochs

clc
clear
format long g

% Julian Date.
d = 2457073.05631;
fprintf(1,'%13.5f\n', d);

% Transform into Besselian epoch.
e = iauEpb ( 0.0, d );
fprintf(1,'B%15.10f\n', e);

% Transform back.
[d1, d2] = iauEpb2jd(e);
fprintf(1,'%17.9f\n', d1+d2);

% The same for Julian epoch.
e = iauEpj(0.0, d);
fprintf (1,'J%15.10f\n', e);
[d1, d2] = iauEpj2jd(e);
fprintf (1,'%17.9f\n', d1+d2);
fprintf (1,'\n');

