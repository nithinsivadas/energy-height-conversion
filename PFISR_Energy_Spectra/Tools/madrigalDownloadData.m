%% Program that can download PFISR data files
% Run madrigalkindatfilter.m before this to obtain expFileArrayStore which
% contains all the file names that use the Baker Code 

load '/home/nithin/NithinBU/Aug 3 2015/expFileArratStore08_14.mat'

cgiurl='http://isr.sri.com/madrigal/cgi-bin/';
[d1, d2] = size(expFileArrayStore);

for i=1:1:d2

fileNameStr=sprintf('DataFile_%s_%d',expFileArrayStore(i).name(27:30),i);
result=madDownloadFile(cgiurl, expFileArrayStore(i).name,fileNameStr,'Nithin Sivadas','nithin@bu.edu','Boston University','hdf5');

end;
