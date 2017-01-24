function [dataMSP] = extract_msp_data( fileNameStr, timeMinStr, timeMaxStr)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

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


