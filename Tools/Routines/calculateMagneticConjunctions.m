%% Routine to calculate the times of conjunctions

%% Initialize
clear all;
timeMinStr = '01 Mar 2008 8:00:00';
timeMaxStr = '15 April 2008 12:00:00';
groundstation.GDZ = [65.126, -147.47];
groundstation.field = 'pokerflat';
conjunction.probeStr = 'pokerflat';
conjunction.radius = 200; %km
f= filesep();
% thmProbeStr = 'tha,thb,thc,thd,the';
thmProbeStr = 'the,thd';
stopAlt = 110; % km
hemiFlag = 1; % Northern Hemisphere
% magFieldModel = 7; %Tsyganenko 1996
magFieldModel = 10;
h5FilePath = [initialize_root_path(),'LargeFiles',f,'magConjunctions',f];
h5FileStr  = ([datestr(timeMinStr,'YYYYmmDD'),...
    'to',datestr(timeMaxStr,'YYYYmmDD'),'_',...
    char(strjoin(string(strsplit(thmProbeStr,',')),'_')),'.h5']);
multiWaitbar('CloseAll');
%% Creating a time array with cell values at every month
timeMinVec = datevec(timeMinStr);
timeMaxVec = datevec(timeMaxStr);
timeMin = datenum(timeMinStr);
timeMax = datenum(timeMaxStr);

thisTime = timeMin;
thisTimeVec = timeMinVec;
iMonth = 1;
while thisTime<timeMax
    downloadTimes(iMonth) = thisTime;
    if thisTimeVec(2)<12                    % Incrementing by one month
        thisTimeVec(2) = thisTimeVec(2)+1;
    else
        thisTimeVec(1) = thisTimeVec(1)+1;
        thisTimeVec(2) = 1;
    end
    thisTime = datenum(thisTimeVec);
    iMonth=iMonth+1;
end
downloadTimes(end) = timeMax;

%% Calculating Magnetic conjunctions
nMonth = length(downloadTimes);
[roothPath,comp] = initialize_root_path;
tempStorage = [initialize_root_path,'..',filesep,'Temp',filesep];

minute = 1/(24*60); % 1 min is 1/(24*60) days
multiWaitbar('Calculating magnetic conjunctions',0);

% Creating H5 file to store data
origPath = pwd;
if ~isfolder(h5FilePath)
    mkdir(h5FilePath);
end
cd(h5FilePath);
if isfile([h5FilePath,h5FileStr])
    warning(['Overwriting : ',h5FileStr]);
else
    h5create(h5FileStr,'/conjunctions/flag',[inf inf],'Datatype','single','ChunkSize',[10 50000]);
    h5create(h5FileStr,'/conjunctions/time',[inf inf],'Datatype','double','ChunkSize',[1 50000]);
end
cd(origPath);
h5Pointer = 1;
for iMonth=1:1:nMonth
    % Generating the input data for find_magnetic_conjunctions
    multiWaitbar('Calculating magnetic conjunctions','Increment',1/nMonth,...
        'Relabel',['Calculating mag. conj. for ',datestr(downloadTimes(iMonth),'mmm YYYY')]);
    %Downloads themis and omni date for the month
    omniData = process_omni_data(datestr(downloadTimes(iMonth)),tempStorage);
    themisData = process_themis_data(datestr(downloadTimes(iMonth)),tempStorage,thmProbeStr);
    % Generating minute-wise time array 
    if iMonth ==1
        if nMonth ~=1
            thisTimeArray = downloadTimes(iMonth):minute:downloadTimes(iMonth+1);
        else
            thisTimeArray = timeMin:minute:timeMax;
        end
    elseif iMonth == nMonth
        startVec = timeMaxVec;
        startVec(3)=1;
        startTime = datenum(startVec);
        thisTimeArray = startTime:minute:downloadTimes(iMonth);
    else
        thisTimeArray = omniData.minutely.time';
    end
    % Store the coordinate system and time into data structure usable by
    % find_magnetic_conjunctions
    themisData=rmfield(themisData,'info');
    scField=strsplit(thmProbeStr,',');
    for iSC = 1:1:length(scField)
        thisProbe = char(scField(iSC));
        spacecraft.(thisProbe).GEO = themisData.(thisProbe).state.XYZ_GEO;
        spacecraft.(thisProbe).time = themisData.(thisProbe).state.time;
    end
    % Find magneticfield conjunctions 
    probesOne = find_magnetic_conjunctions(thisTimeArray,omniData,spacecraft,...
        groundstation,conjunction,stopAlt,hemiFlag,magFieldModel);
    probesOne.time = thisTimeArray;
    probes(iMonth) = probesOne;
        
% Listing down the conjunction start and end times 
    multiWaitbar('Recording conjunction start and end times',0,'Color',[0.2 0.9 0.3]);
    for iProbe = 1:1:length(probes(iMonth).probeNames)
        multiWaitbar('Recording conjunction start and end times',...
            'Increment',1/length(probes(iMonth).probeNames));
        thisProbe = char(probes(iMonth).probeNames(iProbe));
        if ~strcmp(thisProbe,conjunction.probeStr)
            startIndx = find(diff(probes(iMonth).flag(:,iProbe))==1);
            stopIndx  = find(diff(probes(iMonth).flag(:,iProbe))==-1);
            % Appending arrays in dataConjunctions after each month
            % iteration
            if iMonth==1
                dataConjunctions.(thisProbe).startTime = [];
                dataConjunctions.(thisProbe).stopTime = [];
                dataConjunctions.(thisProbe).startGDZ = [];
                dataConjunctions.(thisProbe).stopGDZ = [];
            end
            dataConjunctions.(thisProbe).startTime =...
                [dataConjunctions.(thisProbe).startTime;thisTimeArray(startIndx)'];
            dataConjunctions.(thisProbe).stopTime =...
                [dataConjunctions.(thisProbe).stopTime; thisTimeArray(stopIndx)'];
            dataConjunctions.(thisProbe).startGDZ =...
                [dataConjunctions.(thisProbe).startGDZ; probes(iMonth).(thisProbe).GDZ(startIndx,:)];
            dataConjunctions.(thisProbe).stopGDZ =... 
                [dataConjunctions.(thisProbe).stopGDZ; probes(iMonth).(thisProbe).GDZ(stopIndx,:)];
        end
    end
    
    multiWaitbar('Recording conjunction start and end times','Reset');
    multiWaitbar(['Calculating mag. conj. for ',datestr(downloadTimes(iMonth),'mmm YYYY')],...
        'Relabel','Calculating magnetic conjunctions');
    dataSize=size(probes(iMonth).flag');
    timeSize=dataSize(2);
    h5write([h5FilePath,h5FileStr],'/conjunctions/flag',single(probes(iMonth).flag)',[1 h5Pointer],dataSize);
    h5write([h5FilePath,h5FileStr],'/conjunctions/time',probes(iMonth).time,[1 h5Pointer],[1 timeSize]);
    h5Pointer = h5Pointer + dataSize(2);
end   
