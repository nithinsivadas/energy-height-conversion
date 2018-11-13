function [p] = combine_2D_plots_v2(inputH5FileStr,figureHandle,...
    varargin)
%combine_2D_plots Plot optical images, energy flux maps, and or magnetic
%field maps over one another.

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
expectedMaps = {'OpticalImage','EnergyFluxMap','MagneticFieldMap','NoMap'};
expectedSites = {'gako','fykn','mcgr','whit','inuv','kian','dasc','pokerFlat'};
expectedMagMapContour = {'RE','Lm','Lstar','Kc'}; 
expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};

% Map type, should be a cell array
addParameter(p,'maps',{'OpticalImage','EnergyFluxMap'},@(x) iscell(x));

% Corresponding instrument-site or magnetic field model cell array
addParameter(p,'sites',{'pokerFlat','pokerFlat'},@(x) iscell(x));

addParameter(p,'plotContours','RE',@(x) any(validatestring(x,expectedMagMapContour)));
addParameter(p,'thisTime',datenum('26 Mar 2008 11:00'),validScalarPosNum);
addParameter(p,'magneticFieldModel','TS96',@(x) any(strcmp(x,expectedMagFieldModels)));

addParameter(p,'peakIonizationAltitude',nan, validScalarPosNum); %Reduce initial processing time of time independent coordinates of energyFlux maps
addParameter(p,'showComments',false,@(x) islogical(x));

addParameter(p,'latLim',[63 67]);
addParameter(p,'lonLim',[-153 -143]);
addParameter(p,'elCutoff',30,validScalarPosNum);
addParameter(p,'transparency',0.8,validScalarPosNum);
addParameter(p,'setOpticalLabel',false, @(x) islogical(x));
addParameter(p,'imageSize', 512, validScalarPosNum);

addParameter(p, 'opticalLim', [300 600]);
addParameter(p, 'deltaLat', 1);
addParameter(p, 'deltaLon', 5);

addParameter(p, 'eFluxLim', [8 10]);

addParameter(p,'energySlice',100,validScalarPosNum);
addParameter(p,'setEnergyTimeLabelOn',true);

addParameter(p,'contourLineArray',1:1:30);
addParameter(p,'contourLabelArray',[1,5:5:30]);
addParameter(p,'setMagneticFieldTimeLabelOn',false);
addParameter(p,'setFieldLabelOn',false);

addParameter(p,'figureLength',148,validScalarPosNum);
addParameter(p,'figureBreadth',210,validScalarPosNum);

addParameter(p, 'imageStoreDir',['.',filesep,'TemporaryImages',filesep]);
addParameter(p, 'setStoreImage',true);

addRequired(p,'inputH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));
addRequired(p,'figureHandle',@(x)isfigure(x));


parse(p,inputH5FileStr,figureHandle,varargin{:});

nTime = length(p.Results.thisTime);

% We have three types of maps: 1. Optical, 2. EnergyFlux, 3. MagneticFieldLines
% Optical - layer 1, EnegryFlux - layer 2, MagneticField - layer 3
maps.ID = 1:1:length(p.Results.maps);
maps.optical = find(strcmp(p.Results.maps,expectedMaps{1}));
maps.energy = find(strcmp(p.Results.maps,expectedMaps{2}));
maps.magneticField = find(strcmp(p.Results.maps,expectedMaps{3}));

nOptical = length(maps.optical);
nEnergy = length(maps.energy);
nMagneticField = length(maps.magneticField);

% Loading non-time-dependent variables
comment('Loading time independent data...',p.Results.showComments);
if nOptical>0
    for i=1:1:nOptical
        maps.opticalData(i) =  get_2D_plot_inputs_time_independent(inputH5FileStr,...
        'plotModeStr','OpticalImage','site',p.Results.sites{maps.optical(i)});
    end
end
if nEnergy>0
    for i=1:1:nEnergy
        maps.energyData(i) =  get_2D_plot_inputs_time_independent(inputH5FileStr,...
        'plotModeStr','EnergyFluxMap','site',p.Results.sites{maps.energy(i)},...
        'peakIonizationAltitude',p.Results.peakIonizationAltitude);
    end
end
if nMagneticField>0
    for i=1:1:nMagneticField
        maps.magneticFieldData(i) =  get_2D_plot_inputs_time_independent(inputH5FileStr,...
        'plotModeStr','MagneticFieldMap','site',p.Results.sites{maps.magneticFieldData(i)});
    end
end

if  p.Results.setStoreImage == true
    if ~isfolder(p.Results.imageStoreDir)
        comment('Creating image storage directory...',p.Results.showComments);
        mkdir(p.Results.imageStoreDir);
    end
end

for iTime = 1:1:nTime
    imageName = strcat('figure_',datestr(p.Results.thisTime(iTime),'HH_MM_SS'));
    if isfile(strcat(p.Results.imageStoreDir,imageName,'.png')) && p.Results.setStoreImage == true
    warning('Image file already exists, skipping calculation...');
        % Create two axes
    else
        resize_figure(figureHandle,p.Results.figureLength,p.Results.figureBreadth); %A5 Paper Size, 148 cm vert, 210 cm horizontal
        
        if nOptical>0 %Plot opitcal images on a map
            comment('Plotting optical images...',p.Results.showComments);
            for i=1:1:nOptical
            thisTimeIndx = find_time(maps.opticalData(i).time,datestr(p.Results.thisTime(iTime)));
            dascData = get_2D_plot_inputs_at_time(inputH5FileStr,...
                    'plotModeStr',expectedMaps{1},...
                    'plotData',maps.opticalData(i),...
                    'thisTimeIndx', thisTimeIndx,...
                    'site', p.Results.sites{maps.optical(i)});
            
            if strcmpi(p.Results.sites(maps.optical(i)),'pokerFlat')
                intensityScale = 250./(dascData.background-25); % Custom scaling
            else
                intensityScale = 250./dascData.background; %Scale size of the intensities measured
            end
            % Cutting-off elevation below the specified
            dascData.image(maps.opticalData(i).elevation<p.Results.elCutoff) = nan;
            dascData.image(isnan(maps.opticalData(i).elevation))=nan;
            if sum(~isnan(dascData.image(:)))>0
%             if strcmp(p.Results.sites{maps.optical(i)},'pokerFlat')
%                 dascData.image((dascData.image(:)>=30000))=nan;
%             end
                if i==1
                dascData.image(dascData.image>=65536)=nan; % Remove saturated pixels (16Bit)
                %Correcting all cameras (rudimentary method); nightsky background intensity
                
%                 maxMedianMarker = nightsky(i)./max(dascData.image(:));
                axesHandleOptical=axes;
                [axesmHandleOptical, hOptical(i)] = plot_DASC_geodetic((dascData.image(:)').*intensityScale,...
                    dascData.thisTime, dascData.latitude(:)', dascData.longitude(:)',...
                    p.Results.imageSize, p.Results.latLim, p.Results.lonLim, ...
                    p.Results.deltaLat,p.Results.deltaLon);
                    colormap(axesHandleOptical,'viridis');
                    cbOptical = colorbar(axesHandleOptical,'eastoutside');
                    ylabel(cbOptical,'[a.u.]');
                    caxis(axesHandleOptical,p.Results.opticalLim);
                    
%                     hold on;
%                     plotm([60, 59, 58],[-150, -147, -145],'*');
                else
                hold on;
%                 dascData.image(dascData.image>=65536)=nan; % Remove saturated pixels (16Bit)
                
                [~,hOptical(i)]=plot_DASC_geodetic((dascData.image(:)').*intensityScale,...
                    dascData.thisTime, dascData.latitude(:)', dascData.longitude(:)',...
                    p.Results.imageSize, p.Results.latLim, p.Results.lonLim, ...
                    p.Results.deltaLat,p.Results.deltaLon);
                    colormap(axesHandleOptical,'viridis');
                    caxis(axesHandleOptical,p.Results.opticalLim);
                
                end
                if p.Results.setOpticalLabel
                    textm(maps.opticalData(i).sensorloc(1), maps.opticalData(i).sensorloc(2),...
                       {char(upper(p.Results.sites(maps.optical(i)))),[datestr(dascData.thisTime,'HH:MM:SS'),' UT']});
                end
                alpha(hOptical(i),p.Results.transparency);         
            else
                if i==1
                axesHandleOptical=axes;
                end
            end
            end
        end
        
        if nEnergy>0 % Plot energyFluxMaps
            comment('Plotting energy flux maps...',p.Results.showComments);
            for i=1:1:nEnergy
                thisTimeIndx = find_time(maps.energyData(i).time,datestr(p.Results.thisTime(iTime)));
                
                pfisrData = get_2D_plot_inputs_at_time(inputH5FileStr,...
                    'plotModeStr',expectedMaps{2},...
                    'plotData',maps.energyData(i),...
                    'thisTimeIndx', thisTimeIndx);
                
                
                latWidth = p.Results.latLim(2)-p.Results.latLim(1);
                lonWidth = p.Results.lonLim(2)-p.Results.lonLim(1);
                pfisrData.diffEnergyFlux(imag(pfisrData.diffEnergyFlux(:))~=0)=nan;
                
                if i==1
                axesHandleEnergy=axes;
                [axesmHandleEnergy,hEnergy]=plot_2D_energy_slice_geodetic_v2018...
                    (pfisrData.diffEnergyFlux,...
                pfisrData.latitude, pfisrData.longitude,...
                pfisrData.zEnergyBin, pfisrData.thisTime,...
                p.Results.energySlice,latWidth,lonWidth,true,...
                p.Results.setEnergyTimeLabelOn);  
                colormap(axesHandleEnergy,'inferno');
                cbEnergy = colorbar(axesHandleEnergy,'westoutside');
                ylabel(cbEnergy,'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
                caxis(axesHandleEnergy,p.Results.eFluxLim);
%                 hold on;
%                 plotm([60, 59, 58],[-150, -147, -145],'o');
                else
                hold on;
                plot_2D_energy_slice_geodetic_v2018...
                    (pfisrData.diffEnergyFlux,...
                pfisrData.latitude, pfisrData.longitude,...
                pfisrData.zEnergyBin, pfisrData.thisTime,...
                p.Results.energySlice,latWidth,lonWidth,true,...
                p.Results.setEnergyTimeLabelOn);
                end             
            end
        end

        if nMagneticField>0
            comment('Plotting magnetic field footprints...',p.Results.showComments);
            for i = 1:1:nMagneticField
            thisTimeIndx = find_time(maps.magneticFieldData(i).time,datestr(p.Results.thisTime(iTime)));
            magFieldData = get_2D_plot_inputs_at_time(inputH5FileStr,...
                    'plotModeStr',expectedMaps{3},...
                    'plotData',maps.magneticFieldData(i),...
                    'thisTimeIndx', thisTimeIndx,...
                    'magFieldModelStr',p.Results.magneticFieldModel,...
                    'energySlice',p.Results.energySlice,...
                    'getKc',strcmp(p.Results.plotContours,'Kc'));
                
                    if i==1
                    axesHandleMagneticField=axes;
                    
                    if strcmp(p.Results.plotContours,'RE')
                         plotVariable = magFieldData.RE;
                    elseif strcmp(p.Results.plotContours,'Lm')
                        plotVariable = magFieldData.Lm;
                    elseif strcmp(p.Results.plotContours,'Lstar')
                        plotVariable = magFieldData.Lstar;
                    elseif strcmp(p.Results.plotContours,'Kc')
                        plotVariable = abs(magFieldData.Kc);
                    end

                    [axesmHandleMagneticField, hMagnetic, cMagnetic] = plot_2D_magnetic_foot_points...
                        (plotVariable, magFieldData.ionosphereCoord,...
                        'plotVariableName',p.Results.plotContours,...
                        'BfieldModelStr',magFieldData.magFieldModelStr, 'thisTimeBfieldModel',magFieldData.thisTime,...
                        'setMapOn',true, 'latLim', p.Results.latLim,'lonLim', p.Results.lonLim,'contourLineArray', p.Results.contourLineArray,...
                        'setFieldLabelOn',false,'setTimeLabelOn', p.Results.setMagneticFieldTimeLabelOn);
                    else
                        hold on;
                    
                    if strcmp(p.Results.plotContours,'RE')
                         plotVariable = magFieldData.RE;
                    elseif strcmp(p.Results.plotContours,'Lm')
                        plotVariable = magFieldData.Lm;
                    elseif strcmp(p.Results.plotContours,'Lstar')
                        plotVariable = magFieldData.Lstar;
                    elseif strcmp(p.Results.plotContours,'Kc')
                        plotVariable = abs(magFieldData.Kc);
                    end

                    plot_2D_magnetic_foot_points...
                        (plotVariable, magFieldData.ionosphereCoord,...
                        'plotVariableName',p.Results.plotContours,...
                        'BfieldModelStr',magFieldData.magFieldModelStr, 'thisTimeBfieldModel',magFieldData.thisTime,...
                        'setMapOn',true, 'latLim', p.Results.latLim,'lonLim', p.Results.lonLim,'contourLineArray', p.Results.contourLineArray,...
                        'setFieldLabelOn',false,'setTimeLabelOn', p.Results.setMagneticFieldTimeLabelOn);
                    end 
            end     
        end
        
        comment('Linking Axes...',p.Results.showComments);
        if nOptical > 0 && nEnergy > 0
            linkaxes([axesHandleOptical,axesHandleEnergy]);
            set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'),'XLim',get(axesHandleOptical,'XLim'));
            if nMagneticField > 0
                tempAxesHandle.XLim = axesHandleMagneticField.XLim;
            end
    

            setm(axesmHandleEnergy,'MapProjection','lambertstd',...
                'MapLatLimit',getm(axesmHandleOptical,'MapLatLimit'),...
                'MapLonLimit',getm(axesmHandleOptical,'MapLonLimit'),...
            'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
            'PLineVisible','off','MLineVisible','off');  
            setm(axesHandleEnergy,'Origin',getm(axesHandleOptical,'Origin'),...
                'FLatLimit',getm(axesmHandleOptical,'FLatLimit'),...
                'FLonLimit',getm(axesmHandleOptical,'FLonLimit'));
            alpha(axesHandleEnergy, 0.5);
        
        elseif nMagneticField > 0 && nEnergy > 0
            linkaxes([axesHandleMagneticField,axesHandleEnergy]);
            set(axesHandleEnergy,'Position',get(axesHandleMagneticField,'Position'),'XLim',get(axesHandleMagneticField,'XLim'));
            
            tempAxesHandle.XLim = axesHandleMagneticField.XLim;
            
            setm(axesmHandleEnergy,'MapProjection','lambertstd',...
                'MapLatLimit',getm(axesmHandleMagneticField,'MapLatLimit'),...
                'MapLonLimit',getm(axesmHandleMagneticField,'MapLonLimit'),...
            'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
            'PLineVisible','off','MLineVisible','off');        
            setm(axesHandleEnergy,'Origin',getm(axesHandleMagneticField,'Origin'),...
                'FLatLimit',getm(axesmHandleMagneticField,'FLatLimit'),...
                'FLonLimit',getm(axesmHandleMagneticField,'FLonLimit'));
            alpha(axesHandleEnergy, 0.5);    
           
        elseif nOptical > 0 && nEnergy > 0 && nMagneticField > 0
            linkaxes([axesHandleOptical,axesHandleEnergy,axesHandleMagneticField]);
            set(axesHandleEnergy,'Position',get(axesHandleOptical,'Position'),'XLim',get(axesHandleOptical,'XLim'));
            set(axesHandleMagneticField,'Position',get(axesHandleOptical,'Position'),'XLim',get(axesHandleOptical,'XLim'));

            tempAxesHandle.XLim = axesHandleMagneticField.XLim;
            
            setm(axesmHandleEnergy,'MapProjection','lambertstd','MapLatLimit',getm(axesmHandleOptical,'MapLatLimit'),...
                'MapLonLimit',getm(axesmHandleOptical,'MapLonLimit'),...
            'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
            'PLineVisible','off','MLineVisible','off');

            setm(axesmHandleMagneticField,'MapProjection','lambertstd','MapLatLimit',getm(axesmHandleOptical,'MapLatLimit'),...
                'MapLonLimit',getm(axesmHandleOptical,'MapLonLimit'),...
            'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
            'PLineVisible','off','MLineVisible','off');
            setm(axesmHandleEnergy,'Origin',getm(axesHandleOptical,'Origin'));
            setm(axesHandleEnergy,'FLatLimit',getm(axesmHandleOptical,'FLatLimit'),...
                'FLonLimit',getm(axesmHandleOptical,'FLonLimit'));
            setm(axesHandleMagneticField,'Origin',getm(axesmHandleEnergy,'Origin'));
            setm(axesHandleMagneticField,'FLatLimit',getm(axesmHandleOptical,'FLatLimit'),...
                'FLonLimit',getm(axesmHandleOptical,'FLonLimit'));

            alpha(axesHandleEnergy, 0.5);
            alpha(axesHandleMagneticField, 0.5);
        end
        
        comment('Making some adjustments on final figure...',p.Results.showComments);
        if nMagneticField > 0 && p.Results.setFieldLabelOn
        htext = clabelm(cMagnetic,hMagnetic,p.Results.contourLabelArray,...
            'LabelSpacing',723);
        set(htext,'BackgroundColor','cyan','margin',2);
        end

        if nMagneticField > 0
            set(axesHandleMagneticField,'XLim',tempAxesHandle.XLim);
        end
        
        if p.Results.setStoreImage == true
            comment('Storing the image...',p.Results.showComments);
            export_fig(strcat(p.Results.imageStoreDir,imageName,'.png'),'-r300','-png','-nocrop');
            close(figureHandle);
        end
    
    end
title(datestr(p.Results.thisTime(iTime)));
comment('Done...',p.Results.showComments);       
end
   
end

function comment(inputStr,showComments)
    if showComments || nargin<2
        disp(inputStr);
    end
end

function OK = isfigure(h)
 if strcmp(get(h,'type'), 'figure')
     OK = true;
 else
     OK = false;
 end
     
end

% function minElFilter = elevation_cut_off(dascData,minEl)
%     if nargin<2
%         minEl = 30;
%     end
%     minElFilter = ones(size(dascData.elevation));
%     minElFilter(find(dascData.elevation<minEl)) = 0;
% end