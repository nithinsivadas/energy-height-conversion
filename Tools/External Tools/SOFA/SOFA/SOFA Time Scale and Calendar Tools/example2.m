% example 2: Formatting conventions

clc
clear
format long g

% The time.
ih = 23;
im = 5;
s = 11.630799;
fprintf (1,'%2d:%2.2d:%9.6f\n', ih, im, s);

% Express as a fraction of 1 day.
f = iauTf2d ('+', ih, im, s);
fprintf(1,'%14.12f\n', f );

% Six hours earlier.
f = f - 0.25;

% Report to 1 ms precision.
[~, ihmsf] = iauD2tf(3, f);
fprintf(1,'%2d:%2.2d:%2.2d.%3.3d\n', ihmsf(1),ihmsf(2),ihmsf(3),ihmsf(4));
fprintf(1,'\n');

