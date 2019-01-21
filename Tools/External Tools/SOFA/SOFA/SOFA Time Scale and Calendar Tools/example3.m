% example 3: Julian date

clc
clear
format long g

% Date and time.
iy = 2008; im = 2; id = 29;
ihour = 23; imin = 59; sec = 59.9;
fprintf(1,'%4d/%2.2d/%2.2d%3d:%2.2d:%4.1f\n',iy, im, id, ihour, imin, sec);
% Express as 2-part JD.
[d1, d2] = iauCal2jd ( iy, im, id);
d = iauTf2d('+', ihour, imin, sec);
d2 = d2 + d;
fprintf(1,'%9.1f +%13.6f =%15.6f\n', d1, d2, d1 + d2);

% Express as calendar date and fraction of a day.
[iy, im, id, fd] = iauJd2cal(d1, d2);

d = id + fd;
fprintf(1,'%4d/%2.2d/%9.6f\n', iy, im, d);

% Round to 0.001 day.
iymdf = iauJdcalf(3, d1, d2);
fprintf(1,'%4d/%2.2d/%2.2d.%3.3d\n',iymdf(1),iymdf(2),iymdf(3),iymdf(4));
fprintf(1,'\n');

