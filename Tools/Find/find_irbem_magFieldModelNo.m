function magFieldModelNo = find_irbem_magFieldModelNo(magFieldModelStr)
%find_irbem_magFieldModelNo: Finds the magnetic field number corresponding
%to the strings 'NoExternalField','MF75','TS87short','TS87long',...
%     'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
%     'TS04storm','Alexeev2000'

expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};
if sum(strcmpi(expectedMagFieldModels,magFieldModelStr))==0
    error('Incorrect magnetic field model string');
end
magFieldModelNo=find(strcmpi(expectedMagFieldModels,magFieldModelStr))-1;
end

