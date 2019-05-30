MyPath = 'C:\Users\nithin\Documents\GitHub\energy-height-conversion\Tools\';
cd(MyPath);
files = dir('**/*.m');
filesCell = struct2cell(files);

% Find
Convert = filesCell(1,contains(filesCell(2,:),'\Convert'))';
ExternalTools = filesCell(1,contains(filesCell(2,:),'\External Tools'))';
Find = filesCell(1,contains(filesCell(2,:),'\Find'))';
Get = filesCell(1,contains(filesCell(2,:),'\Get'))';
Plot = filesCell(1,contains(filesCell(2,:),'\Plot'))';
Read = filesCell(1,contains(filesCell(2,:),'\Read'))';
Routines = filesCell(1,contains(filesCell(2,:),'\Routines'))';
% MISC
Misc = filesCell(1,contains(filesCell(2,:),'\Misc'))';



