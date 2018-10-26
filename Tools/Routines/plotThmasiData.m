
%% THEMIS ASI plotting script
clear all;

%% Initializing
thmasiCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\ThemisASI';
thmasiCDFName = 'thg_l1_asf_gako_2008032611_v01.cdf';

[data,info] = cdfread([thmasiCDFPath,filesep,thmasiCDFName]);

% thmasiCalFile = 'thg_l1_asf_fykn_00000000_v01.cdf';
thmasiCalFile1 = 'thg_l2_asc_gako_19700101_v01.cdf';
% thmasiCalFile2 = 'themis_skymap_gako_20070401.cdf';

[dataCal,infoCal] = cdfread([thmasiCDFPath,filesep,thmasiCalFile]);
[infoCal2] = spdfcdfinfo([thmasiCDFPath,filesep,thmasiCalFile1]);