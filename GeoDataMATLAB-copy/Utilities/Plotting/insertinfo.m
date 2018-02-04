function stroutall = insertinfo(strin,varargin)
% This function will augment input strings or cell arrays of strings.
% by John Swoboda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramstr = varargin(1:2:end);
paramvals = varargin(2:2:end);
poss_labels={'key','time','timend'};
varnames = {'key','posix','posixend'};
vals = {'',nan,nan};
checkinputs(paramstr,paramvals,poss_labels,vals,varnames);

cellin = false;
if iscell(strin)
    stroutall = cell(size(strin));
    cellin = true;
else
    strin = {strin};
end

for k =1:length(strin);

    strout = strrep(strin{k},'$k',key);
    if isnan(posix);
        strout=strrep(strout,'$tu','');
        strout=strrep(strout,'$tdu','');
    else
        curdt = unixtime2matlab(posix);
        curdte = unixtime2matlab(posixend);
        markers = {
            '$thmsehms',...%UT hours minutes seconds - hours minutes seconds
            '$thmehm',...%UT hours minutes - hours minutes
            '$tmsems',...%UT minutes seconds - minutes seconds
            '$thms',...%UT hours minutes seconds
            '$thm',...%UT hours minutes
            '$tms',...%UT minutes seconds
            '$tmdyhms',...%UT month/day/year hours minutes seconds
            '$tmdyhm',...%UT month/day/year hours minutes
            '$tmdhm'...%UT month/day hours minutes
            };
        datestrcell = {...
            [datestr(curdt,'HH:MM:SS'),' - ',datestr(curdte,'HH:MM:SS'),' UT'],...
            [datestr(curdt,'HH:MM'),' - ',datestr(curdte,'HH:MM'),' UT'],...
            [datestr(curdt,'MM:SS'),' - ',datestr(curdte,'MM:SS'),' UT'],...
            [datestr(curdt,'HH:MM:SS'),' UT'],...
            [datestr(curdt,'HH:MM'),' UT'],...
            [datestr(curdt,'MM:SS'),' UT'],...
            [datestr(curdt,'mm/dd/yyyy HH:MM:SS'),' UT'],...
            [datestr(curdt,'mm/dd/yyyy HH:MM'),' UT'],...
            [datestr(curdt,'mm/dd HH:MM'),' UT'],...
            };
        for imark =1:length(markers)
            strout=strrep(strout,markers{imark},datestrcell{imark});
        end   
    end
    if cellin
        stroutall{k} = strout;
    else
        stroutall = strout;
    end
end

    