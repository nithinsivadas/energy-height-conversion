function [themisData, matFilePath] = process_themis_data(dateStr,dataStoreDir,probeStr,dataType)
%UNTITLED2 Summary of this function goes here
%   Process and download themis data for a month
f=filesep;
 if nargin < 4
        dataType = 'state';
    end
    if nargin < 3
        probeStr = 'tha,thb,thc,thd,the';
    end
    if nargin < 2
        dataStoreDir = [initialize_root_path,'LargeFiles',f];
    end
    if nargin < 1
        dateStr = '26 Mar 2008';
    end
    
    probe = strsplit(probeStr,',');
    date = datenum(dateStr);
    if strcmp(dataType,'state')
    cdfLocalDir = [dataStoreDir,'themis',f,'state',f];
     for thisSC = 1:length(probe)
        thisProbe = char(probe(thisSC));
        cdfFileStr = [thisProbe,'_or_ssc_',datestr(date,'yyyymm'),'01_v01.cdf'];
        if ~isfile([cdfLocalDir,cdfFileStr])
            warning('Downloading CDF from online SPDF database');
            download_themis(dateStr,dataStoreDir,probeStr,dataType);
        end
        
        varNames = [{'XYZ_GSM'};{'XYZ_GEO'}];
        themisData.(thisProbe).state.time = spdfcdfread...
            ([cdfLocalDir,cdfFileStr],'Variables',{'Epoch'},...
            'ConvertEpochToDatenum', true);
        themisData.(thisProbe).state.XYZ_GSM = spdfcdfread(...
            [cdfLocalDir,cdfFileStr],'Variables',{'XYZ_GSM'});
        themisData.(thisProbe).state.XYZ_GEO = spdfcdfread(...
            [cdfLocalDir,cdfFileStr],'Variables',{'XYZ_GEO'});

        info = spdfcdfinfo([cdfLocalDir,cdfFileStr]);
        [tf, loc] = ismember(info.VariableAttributes.UNITS(:,1),varNames');
        [~, p] = sort(loc(tf));
        VarNamesID = find(tf);
        VarNamesID = VarNamesID(p);
        themisData.(thisProbe).state.Units = info.VariableAttributes.UNITS(VarNamesID,:);


        [tf, loc] = ismember(info.VariableAttributes.FIELDNAM(:,1),varNames');
        [~, p] = sort(loc(tf));
        VarNamesID = find(tf);
        VarNamesID = VarNamesID(p);
        themisData.(thisProbe).state.fieldNames = info.VariableAttributes.FIELDNAM(VarNamesID,:);

     end 
    end
    
themisData.info = ['Contains themis data for ',datestr(date,'mmm YYYY')];

%% Appending themis data onto themisData.mat file
% Need to check if this works after adding more than one dataType.
themisLocalMatDir = [dataStoreDir,datestr(date,'yyyymmdd'),f];
matFilePath = [themisLocalMatDir,'themisData.mat'];
if isfile(matFilePath)
    themisDataOld=themisData;
    file = load(matFilePath);
    themisData=file.themisData;
    for thisSC = 1:length(probe)
        thisProbe = char(probe(thisSC));
        if isfield(themisData,thisProbe)
            M = [fieldnames(themisData.(thisProbe))' fieldnames(themisDataOld.(thisProbe))';...
                struct2cell(themisData.(thisProbe))' struct2cell(themisDataOld.(thisProbe))'];
            [tmp, rows] = unique(M(1,:), 'last');
            M=M(:, rows);
            themisData.(thisProbe)=struct(M{:});
        else
            themisData.(thisProbe)=themisDataOld.(thisProbe);
        end
    end
    themisData.info = themisDataOld.info;
    save(matFilePath,'themisData','-append');
else
    if ~isfolder(themisLocalMatDir)
        mkdir(themisLocalMatDir);
    end 
    save(matFilePath,'themisData');
end

%% Creating a File List that records the the directories where the files are stored
processedFileList = [themisLocalMatDir,'processedFileList.mat'];
if isfile(processedFileList)
    file=load(processedFileList);
    matlabFilePath = file.matlabFilePath;
    matlabFilePath.themis=matFilePath;
    save(processedFileList,'matlabFilePath','-append');
else
    matlabFilePath.themis=matFilePath;
    save(processedFileList,'matlabFilePath');
end

end

