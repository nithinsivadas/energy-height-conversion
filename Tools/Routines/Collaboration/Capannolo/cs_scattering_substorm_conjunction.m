%% Find conjunctions between CSS, REPS and storms, substorms, SCMs
clear all;

fileStr = "G:\My Drive\Research\Projects\Collaborations\Colin Forsyth\sophie_output_90.asc";
forsyth.substormT = extract_data(fileStr,2);

fileStrWalachStorm = "G:\My Drive\Research\Projects\Collaborations\Maria Walach\Walach1981to2019_stormlist.txt";
walach.stormT = extract_data(fileStrWalachStorm,3);

fileStrWalachSCM1 = "G:\My Drive\Research\Projects\Collaborations\Maria Walach\Walach_Milan_2015_final_event_lists_SMCs_without_preceeding_substorms.txt";
walach.SCMT1 = extract_data(fileStrWalachSCM1,4);

fileStrWalachSCM2 = "G:\My Drive\Research\Projects\Collaborations\Maria Walach\Walach_Milan_2015_final_event_lists_SMCs.txt";
walach.SCMT2 = extract_data(fileStrWalachSCM2,4);


% CSSs Luisa
fileStrLuisaCSS = "G:\My Drive\Research\Projects\Collaborations\Capannolo Luisa\Capannolo_2021\clear_CSSs_plots_with_indices_2012_to_2020_appended_quantities_T05_with_midnight_crossing.txt";
luisa.CSST = extract_data(fileStrLuisaCSS,5);

% REPs Luisa
fileStrLuisaREP = "G:\My Drive\Research\Projects\Collaborations\Capannolo Luisa\Capannolo_2021\clear_REPs_plots_with_indices_2012_to_2020_appended_quantities_T05_with_midnight_crossing.txt";
luisa.REPT = extract_data(fileStrLuisaREP,5);

luisa.CSST.substorm_phase = interp1(forsyth.substormT.time,double(forsyth.substormT.phase),luisa.CSST.timeI,'previous'); 
luisa.CSST.substorm_MLT = interp1(forsyth.substormT.time,double(forsyth.substormT.MLT),luisa.CSST.timeI,'previous');
luisa.CSST.substorm_MLAT = interp1(forsyth.substormT.time,double(forsyth.substormT.MLAT),luisa.CSST.timeI,'previous');
luisa.REPT.substorm_phase = interp1(forsyth.substormT.time,double(forsyth.substormT.phase),luisa.REPT.timeI,'previous'); 
luisa.REPT.substorm_MLT = interp1(forsyth.substormT.time,double(forsyth.substormT.MLT),luisa.REPT.timeI,'previous'); 
luisa.REPT.substorm_MLAT = interp1(forsyth.substormT.time,double(forsyth.substormT.MLAT),luisa.REPT.timeI,'previous'); 

%%

for i=1:1:height(luisa.CSST)
    if sum(luisa.CSST.timeI(i)> walach.stormT.timeI  & luisa.CSST.timeI(i)<walach.stormT.timeM)
        luisa.CSST.storm_phase(i) = 1;
        luisa.CSST.storm_SymH_Min(i) = walach.stormT.symH_min(luisa.CSST.timeI(i)> walach.stormT.timeI  & luisa.CSST.timeI(i)<walach.stormT.timeM);
        
    elseif sum(luisa.CSST.timeI(i)> walach.stormT.timeM  & luisa.CSST.timeI(i)<walach.stormT.timeR)
        luisa.CSST.storm_phase(i) = 2;
        luisa.CSST.storm_SymH_Min(i) = walach.stormT.symH_min(luisa.CSST.timeI(i)> walach.stormT.timeM  & luisa.CSST.timeI(i)<walach.stormT.timeR);
        
    elseif sum(luisa.CSST.timeI(i)> walach.stormT.timeR  & luisa.CSST.timeI(i)<walach.stormT.timeE)
        luisa.CSST.storm_phase(i) = 3;
        luisa.CSST.storm_SymH_Min(i) = walach.stormT.symH_min(luisa.CSST.timeI(i)> walach.stormT.timeR  & luisa.CSST.timeI(i)<walach.stormT.timeE);
        
    else
        luisa.CSST.storm_phase(i) = nan;
        luisa.CSST.storm_SymH_Min(i) = nan;
    end
   
end

for i=1:1:height(luisa.REPT)
    if sum(luisa.REPT.timeI(i)> walach.stormT.timeI  & luisa.REPT.timeI(i)<walach.stormT.timeM)
        luisa.REPT.storm_phase(i) = 1;
        luisa.REPT.storm_SymH_Min(i) = walach.stormT.symH_min(luisa.REPT.timeI(i)> walach.stormT.timeI  & luisa.REPT.timeI(i)<walach.stormT.timeM);
        
    elseif sum(luisa.REPT.timeI(i)> walach.stormT.timeM  & luisa.REPT.timeI(i)<walach.stormT.timeR)
        luisa.REPT.storm_phase(i) = 2;
        luisa.REPT.storm_SymH_Min(i) = walach.stormT.symH_min(luisa.REPT.timeI(i)> walach.stormT.timeM  & luisa.REPT.timeI(i)<walach.stormT.timeR);
        
    elseif sum(luisa.REPT.timeI(i)> walach.stormT.timeR  & luisa.REPT.timeI(i)<walach.stormT.timeE)
        luisa.REPT.storm_phase(i) = 3;
        luisa.REPT.storm_SymH_Min(i) = walach.stormT.symH_min(luisa.REPT.timeI(i)> walach.stormT.timeR  & luisa.REPT.timeI(i)<walach.stormT.timeE);
        
    else
        luisa.REPT.storm_phase(i) = nan;
        luisa.REPT.storm_SymH_Min(i) = nan;
    end
end


figure; histogram(luisa.CSST.substorm_phase,'Normalization','pdf'); hold on;histogram(luisa.REPT.substorm_phase,'Normalization','pdf'); xlabel('Substorm Phase'); ylabel('pdf'); legend('CSS','REP');
figure; histogram(luisa.CSST.storm_phase,'Normalization','pdf'); hold on;histogram(luisa.REPT.storm_phase,'Normalization','pdf'); xlabel('Storm Phase'); ylabel('pdf'); legend('CSS','REP');


%%

% outputFile1 = 'G:\My Drive\Research\Projects\Collaborations\Capannolo Luisa\Capannolo_2021\CSS';
% outputFile2 = 'G:\My Drive\Research\Projects\Collaborations\Capannolo Luisa\Capannolo_2021\REP';
% outputFile3 = 'G:\My Drive\Research\Projects\Collaborations\Capannolo Luisa\Capannolo_2021\Readme.txt';
% 
% writetable(luisa.CSST(:,{'sc','datetimeI','datetimeE','substorm_phase','substorm_MLT','substorm_MLAT','storm_phase','storm_SymH_Min'}),[outputFile1,'.xls']);
% writetable(luisa.CSST(:,{'sc','datetimeI','datetimeE','substorm_phase','substorm_MLT','substorm_MLAT','storm_phase','storm_SymH_Min'}),[outputFile1,'.txt'],'Delimiter','tab');
% writetable(luisa.REPT(:,{'sc','datetimeI','datetimeE','substorm_phase','substorm_MLT','substorm_MLAT','storm_phase','storm_SymH_Min'}),[outputFile2,'.xls']);
% writetable(luisa.REPT(:,{'sc','datetimeI','datetimeE','substorm_phase','substorm_MLT','substorm_MLAT','storm_phase','storm_SymH_Min'}),[outputFile2,'.txt'],'Delimiter','tab');
% 
% fid = fopen(outputFile3,'w');
% fprintf(fid,'Data Description \n');
% fprintf(fid,'There are two databases used to generate this data set. \n');
% fprintf(fid,'-----------------------------------------------------------\n');
% fprintf(fid,'Substorms\n');
% fprintf(fid,'-----------------------------------------------------------\n');
% fprintf(fid,['For substorms we use Supplementary Material for "A New Technique for Determining Substorm Onsets and Phases from H-Indices of the Electrojet" (SOPHIE) \n Manuscript 2015JA021343\n',...
% 'Forsyth et al.\n',...
% 'SOPHIE identifies substorms by looking for events about a percentile threshold of dSML/dt  (EPT)  for each year. \n',...
% 'This list has EPT=90, which is a high threshold so only picks out the biggest events (other events may be mis-characterised as ‘quiet/growth’). \n',...
% ' Keep an eye on the MLTs at the time of onset. Ordinarily one would expect onsets to be between 22-04 MLT. Anything outside of that is a bit iffy. \n \n']);
% fprintf(fid,'-----------------------------------------------------------\n');
% fprintf(fid,'Storms\n');
% fprintf(fid,'-----------------------------------------------------------\n');
% fprintf(fid,['For storms we use list generated by Maria-Theresia Walach using the same method as the paper\n',...
% 'SuperDARN observations during geomagnetic storms, geomagnetically active times and enhanced solar wind driving\n',...
% 'by M.-T. Walach & A. Grocott, published in JGR, 2019, doi:10.1029/2019JA026816 \n \n',...
% '-----------------------------------------------------------\n',...
% 'sc: spacecraft \n',...
% 'datetimeI: starting time of the event \n',...
% 'datetimeE: ending time of the event \n',...
% 'substorm_phase: number (1 for growth/energy input, 2 for expansion, 3 for recovery) \n',...
% 'substorm_MLT: magnetic local time of the active SML station at the time \n',...
% 'substorm_MLAT: magnetic latitude of the active SML station at the time \n',...
% 'storm_phase: number (1 for inital phase, 2 for main phase, 3 for recovery phase) \n',...
% 'storm_phase: number (1 for inital phase, 2 for main phase, 3 for recovery phase) \n',...
% 'storm_SymH_Min: is the minimum Sym-H index through the entire storm, and indicator of its strength \n']);
% fclose(fid);





function data=load_ascii_files(loadFile, format, headerlines)
if nargin<3
    headerlines = 1;
end
fileID = fopen(loadFile,'r');
data = textscan(fileID, format, 'headerlines', headerlines);
fclose(fileID);
end

function T1 = extract_data(loadFile, ftype)
    
    if ftype == 1
        format ='%4f %2f %2f %2f %2f %5.2f %5.2f ';
        tempData = load_ascii_files(loadFile, format, 69);
        superMag.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        superMag.time = datenum(datestr(superMag.datetime));
        superMag.mlat = tempData{6};
        superMag.mlt = tempData{7};
    end

    if ftype == 2 %Forsyth SOPHIE substorm phases data
        format = '%4f/%2f/%2f-%2f:%2f:%2f %u %u %5.2f %5.2f';
        tempData = load_ascii_files(loadFile, format, 16);

        T.datetime = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},tempData{6});
        T.time = datenum(datestr(T.datetime));
        T.phase = tempData{7};
        T.flag = tempData{8};
        T.MLT = tempData{9};
        T.MLAT = tempData{10};
    end
    
    if ftype == 3 %Walach Storm data
        format = '%4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f \t %4f';
        tempData = load_ascii_files(loadFile, format, 6);

        T.datetimeI = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        T.datetimeM = datetime(tempData{1+5},tempData{2+5},tempData{3+5},tempData{4+5},tempData{5+5},zeros(size(tempData{5+5})));
        T.datetimeR = datetime(tempData{1+5+5},tempData{2+5+5},tempData{3+5+5},tempData{4+5+5},tempData{5+5+5},zeros(size(tempData{5+5+5})));
        T.datetimeE = datetime(tempData{1+5+5+5},tempData{2+5+5+5},tempData{3+5+5+5},tempData{4+5+5+5},tempData{5+5+5+5},zeros(size(tempData{5+5+5+5})));
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeM = datenum(datestr(T.datetimeM));
        T.timeR = datenum(datestr(T.datetimeR));
        T.timeE = datenum(datestr(T.datetimeE));
        T.symH_min = tempData{21};
    end
    
     if ftype == 4 %Walach SCM data with preceeding pre without preceeding substorms
        format = '%4f%2f%2f %2f:%2f \t %4f%2f%2f %2f:%2f';
        tempData = load_ascii_files(loadFile, format, 15);
        T.datetimeI = datetime(tempData{1},tempData{2},tempData{3},tempData{4},tempData{5},zeros(size(tempData{5})));
        T.datetimeE = datetime(tempData{1+5},tempData{2+5},tempData{3+5},tempData{4+5},tempData{5+5},zeros(size(tempData{5+5})));
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeE = datenum(datestr(T.datetimeE));
        
     end
    
     if ftype == 5 %Luisa, current sheet scattering
        format = ['%s %4f-%2f-%2f/%2f:%2f:%2f.%3f \t %4f-%2f-%2f/%2f:%2f:%2f.%3f',repmat(' %*s',1,44)];
        tempData = load_ascii_files(loadFile, format, 0);
        T.sc = tempData{1};
        T.datetimeI = datetime(tempData{2},tempData{3},tempData{4},tempData{5},tempData{6},tempData{7},tempData{8});
        T.datetimeE = datetime(tempData{1+8},tempData{2+8},tempData{3+8},tempData{4+8},tempData{5+8},tempData{6+8},tempData{7+8});
        T.timeI = datenum(datestr(T.datetimeI));
        T.timeE = datenum(datestr(T.datetimeE));
        
    end
    
    T1 = struct2table(T);
end