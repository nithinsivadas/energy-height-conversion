function [dataMSP] = extract_msp_data( fileNameStr, timeMinStr, timeMaxStr)
%% extract_msp_data Extract Meridional Scanning Photometer Data given Josh's
%msp_vvel.mat file
%-------------------------------------------------------------------------
% Input:
%---------------
% fileNameStr : String containing the file name (Paper22/Data/msp_vvel.mat)
% timeMinStr  : String containing starting time (e.g. '26 Mar 2008 09:00')
% timeMaxStr  : String containing end time (e.g. '26 Mar 2008 11:00')
%-------------------------------------------------------------------------
% Output:
%--------
% dataMSP.intensity5577 : Contains the intenity of 5577A wavelength
% dataMSP.intensity4278 : Contains the intenity of 4278A wavelength
% dataMSP.intensity4861 : Contains the intenity of 4861A wavelength
% dataMSP.intensity6300 : Contains the intenity of 6300A wavelength
% dataMSP.el            : Elevation (90 to -90 deg) y-Axis
% dataMSP.time          : time (in matlab units) x-Axis
% dataMSP.lambda        : Contains the wavelengths of the measured
%                         intensities

%%
%----------------------------------------------------------------------------
% Modified: 7th Feb 2017 
% Created : 21st Dec 2016
% Author  : Nithin Sivadas
% Ref     :
%----------------------------------------------------------------------------
%%

load(fileNameStr);

timeMinNo=find_time(msptime, datenum(timeMinStr));
timeMaxNo=find_time(msptime, datenum(timeMaxStr));

dataMSP.intensity5577=squeeze(diff(:,1,timeMinNo:timeMaxNo));
dataMSP.intensity4278=squeeze(diff(:,2,timeMinNo:timeMaxNo));
dataMSP.intensity4861=squeeze(diff(:,3,timeMinNo:timeMaxNo));
dataMSP.intensity6300=squeeze(diff(:,4,timeMinNo:timeMaxNo));
dataMSP.el = linspace(90,-90, length(diff(:,1,1)));
dataMSP.time = msptime(timeMinNo:timeMaxNo);
dataMSP.lambda = [5577, 4278, 4861, 6300];

dataMSP.intensity5577(dataMSP.intensity5577<=0)=min(dataMSP.intensity5577(dataMSP.intensity5577>0));
dataMSP.intensity4278(dataMSP.intensity4278<=0)=min(dataMSP.intensity4278(dataMSP.intensity4278>0));
dataMSP.intensity4861(dataMSP.intensity4861<=0)=min(dataMSP.intensity4861(dataMSP.intensity4861>0));
dataMSP.intensity6300(dataMSP.intensity6300<=0)=min(dataMSP.intensity6300(dataMSP.intensity6300>0));
end


