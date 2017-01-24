function [ timeAvgData ] = get_time_avg_time_series_data( data, time, timeMin, timeMax )
%GET_TIME_AVG_FLUX Calculates time average flux given the time limits
% Input
% data: Any type of 2D time series data, exmaple: energyFlux [nExnT]
% time: Time array [nTx1]
% timeMin: Minimum time value in matlab units
% timeMax: Maximum time value in matlab units

% Output
% timeAvgData: Time averaged data [nEx1]

%%
%----------------------------------------------------------------------------
% Modified: 26th Sep 2016 
% Created : 26th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%

if(timeMin<time(1) || timeMax>time(end))
    error('Min and Max time values are outside the scope of the time array');
end;

timeMinNo = find_time(time, timeMin);
timeMaxNo = find_time(time, timeMax);

timeAvgData = mean(data(:,timeMinNo:1:timeMaxNo),2);

end

