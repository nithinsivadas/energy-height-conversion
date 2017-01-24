% plotPitchAngleDistributionAlongFieldLine
clear all;

pitchAngle = 0:1:180;
equatorialPitchAngleDistribution = (sin(deg2rad(pitchAngle))).^2;
% equatorialPitchAngleDistribution = ones(1,length(pitchAngle));

L_Shell = 6;

latitude = [0, 20, 40];

for iLatitude = 1:1:length(latitude)
    [j_eq(iLatitude,:), pitchAngleAtLatitude(iLatitude,:), invariantLatitude(iLatitude)] = get_pitch_angle_distribution...
        (equatorialPitchAngleDistribution, pitchAngle,...
        L_Shell, latitude(iLatitude));
end;
iLatitude=iLatitude+1;
latitude(iLatitude) = invariantLatitude(1);
[j_eq(iLatitude,:), pitchAngleAtLatitude(iLatitude,:), invariantLatitude(iLatitude)] = get_pitch_angle_distribution...
        (equatorialPitchAngleDistribution, pitchAngle,...
        L_Shell, latitude(iLatitude));


hFig=figure(1);
resize_figure(hFig);
clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 60; %in mm
p.pack({{panelSize} {panelSize}},{{80}});
p.marginleft = 20;

p(1,1).select();
plot(pitchAngle, equatorialPitchAngleDistribution,'-k');
xlabel('Pitch Angle \alpha_e_q [deg]');
ylabel({'Equatorial Pitch Angle Distribution','[m^-^1sr^-^1s^-^1eV^-^1]]'});
grid on;

p(2,1).select();
iLatitude =1;
plot(pitchAngleAtLatitude(iLatitude,:), j_eq(iLatitude,:),'-k');
hold on;
iLatitude =2;
plot(pitchAngleAtLatitude(iLatitude,:), j_eq(iLatitude,:),'-r');
hold on;
iLatitude =3;
plot(pitchAngleAtLatitude(iLatitude,:), j_eq(iLatitude,:),'-m');
hold on;
iLatitude =4;
plot(pitchAngleAtLatitude(iLatitude,:), j_eq(iLatitude,:),'-g');
grid on;

xlabel('Pitch Angle \alpha_e_q [deg]');
ylabel({'Pitch Angle Distribution','[m^-^1sr^-^1s^-^1eV^-^1]'});
legend(['\lambda = ',num2str(latitude(1)),'^0'],...
    ['\lambda = ',num2str(latitude(2)),'^0'],...
    ['\lambda = ',num2str(latitude(3)),'^0'],...
    ['Loss cone \lambda = ',num2str(latitude(4)),'^0'],...
    'Location','southoutside');
p.pack({{panelSize} {panelSize}},{{80}});
