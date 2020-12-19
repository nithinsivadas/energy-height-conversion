
fileStr = 'G:\My Drive\Research\Projects\Data\20201217-07-37-supermag.txt';
content = fileread(fileStr);


[data,info] = parse_superMag_data(fileStr);

%%

h5OmniStr = 'G:\My Drive\Research\Projects\Data\omni_test.h5';
omni.time = h5read(h5OmniStr,'/Time');
F.SML = griddedInterpolant(data.time,data.SML,'nearest','none');
F.SMU = griddedInterpolant(data.time,data.SMU,'nearest','none');
SML = F.SML(omni.time);
SMU = F.SMU(omni.time);
write_h5_dataset(h5OmniStr,'/Indices/SML',SML',1);
write_h5_dataset(h5OmniStr,'/Indices/SMU',SMU',1);
write_h5_dataset_attribute(h5OmniStr,'/Indices/SML','SML - SuperMag Data Revision 5','[nTimex1]','[nT]');
write_h5_dataset_attribute(h5OmniStr,'/Indices/SMU','SMU - SuperMag Data Revision 5','[nTimex1]','[nT]');
