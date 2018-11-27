%% Add thg to HDF5

%% Initialization
outputFileStr = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-energyFlux_v85.h5';
% sites = {'gako'};
sites={'fykn','inuv','whit','mcgr','kian'};
multiWaitbar('Processing different cameras',0);
for i=1:1:length(sites)
multiWaitbar('Processing different cameras','Increment',1./length(sites));
% outputFileStr ='/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-energyFlux_v85.h5';
siteName = sites{i};
disp(['Processing ',upper(siteName)]);
status = create_thg_hdf5(siteName,outputFileStr);
end
multiWaitbar('Close All');

% %%
% outputFileStr = '/media/nithin/PFISR_002_006/PFISR Processed/Event_List/20080326.001_bc_15sec-energyFlux_v85.h5';
% % sites = {'gako'};
% sites={'fykn','inuv','whit','mcgr','kian'};
% sensorLocationFile = '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/External Tools/thmasi/THEMIS_ASI_Station_List_Nov_2011.xls';
% thgSites = parse_thg_location_xls(sensorLocationFile);
% for i=1:1:length(sites)
%     indx = find(strcmpi(thgSites.code,sites(i)));
%     glat = thgSites.glat(indx);
%     glon = thgSites.glon(indx);
%     input.sensorloc = [glat,glon,0];
%     groupName = ['/',char(upper(thgSites.code(indx))),'/'];
%     dset = input.sensorloc;
%     dset_details.Name = 'sensorloc';
%     dset_details.Location = groupName;
%     attr = '3x1 [lat,lon,alt-meters]';
%     attr_details.Name = 'Dimensions';
%     attr_details.AttachedTo = [groupName,dset_details.Name];
%     attr_details.AttachType = 'dataset';
%     hdf5write(outputFileStr,dset_details,dset,...
%         attr_details,attr,'WriteMode','append');
% end