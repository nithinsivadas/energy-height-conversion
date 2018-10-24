%% Testing phase

poesCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17';
poesCDFName = 'poes_n17_20080326.cdf';

[data,info] = cdfread([poesCDFPath,filesep,poesCDFName]);

poesData.time = cellfun(@(x) todatenum(x),data(:,1));