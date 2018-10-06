function [p] = combine_2D_plots(inputH5FileStr,figureHandle,...
    varargin)
%combine_2D_plots Plot optical images, energy flux maps, and or magnetic
%field maps over one another.

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
expectedMaps = {'OpticalImage','EnergyFluxMap','MagneticFieldMap','NoMap'};
expectedMagMapContour = {'RE','Lm','Lstar'}; 
expectedMagFieldModels = {'NoExternalField','MF75','TS87short','TS87long',...
    'TS89','OP77quiet','OP88dynamic','TS96','OM97','TS01','TS01storm',...
    'TS04storm','Alexeev2000'};
addParameter(p,'map1','NoMap',@(x) any(validatestring(x,expectedMaps)));
addParameter(p,'map2','NoMap',@(x) any(validatestring(x,expectedMaps)));
addParameter(p,'map3','NoMap',@(x) any(validatestring(x,expectedMaps)));

addParameter(p,'map1Data',struct());
addParameter(p,'map2Data',struct());
addParameter(p,'map3Data',struct());

addParameter(p,'plotContours','RE',@(x) any(validatestring(x,expectedMagMapContour)));

addParameter(p,'thisTime',datenum('26 Mar 2008 11:00'),validScalarPosNum);

addParameter(p,'magneticFieldModel','TS96',@(x) any(validatestring(x,expectedMagFieldModels)));

addParameter(p,'latLim',[63 67]);
addParameter(p,'lonLim',[-153 -143]);
addParameter(p,'imageSize', 512, validScalarPosNum);

addParameter(p, 'opticalLim', [300 600]);
addParameter(p, 'eFluxLim', [8 10]);

addParameter(p,'energySlice',100,validScalarPosNum);
addParameter(p,'setEnergyTimeLabelOn',true);

addParameter(p,'contourLineArray',1:1:30);
addParameter(p,'contourLabelArray',[1,5:5:30]);
addParameter(p,'setMagneticFieldTimeLabelOn',true);
addParameter(p,'setFieldLabelOn',true);

addParameter(p,'figureLength',148,validScalarPosNum);
addParameter(p,'figureBreadth',210,validScalarPosNum);

addRequired(p,'inputH5FileStr',@(x)contains(x,{'.h5','.hdf5'}));
addRequired(p,'figureHandle',@(x)isfigure(x));


parse(p,inputH5FileStr,figureHandle,varargin{:});

layer(1)=find(strcmp(p.Results.map1,expectedMaps));
if layer(1)~=4
    thisTimeIndx(1) = find_time(p.Results.map1Data.time,datestr(p.Results.thisTime));
end
layer(2)=find(strcmp(p.Results.map2,expectedMaps));
if layer(2)~=4
    thisTimeIndx(2) = find_time(p.Results.map2Data.time,datestr(p.Results.thisTime));
end
layer(3)=find(strcmp(p.Results.map3,expectedMaps));
if layer(3)~=4
    thisTimeIndx(3) = find_time(p.Results.map3Data.time,datestr(p.Results.thisTime));
end

resize_figure(figureHandle,p.Results.figureLength,p.Results.figureBreadth); %A5 Paper Size, 148 cm vert, 210 cm horizontal

for thisLayer = 1:3
    switch layer(thisLayer)
        case 1 %Plot optical image
            dascData = get_2D_plot_inputs_at_time(inputH5FileStr,...
                'plotModeStr',expectedMaps{1},...
                'plotData',select_map_data(thisLayer,p.Results.map1Data,p.Results.map2Data,p.Results.map3Data),...
                'thisTimeIndx', thisTimeIndx(thisLayer));
            axesHandle(thisLayer)=axes;
            [axesmHandle(thisLayer), hOptical] = plot_DASC_geodetic(dascData.image',...
                dascData.thisTime, dascData.latitude', dascData.longitude',...
                p.Results.imageSize, p.Results.latLim, p.Results.lonLim);
                colormap(axesHandle(thisLayer),'viridis');
                cb(thisLayer) = colorbar(axesHandle(thisLayer),'eastoutside');
                ylabel(cb(thisLayer),'[a.u.]');
                caxis(axesHandle(thisLayer),p.Results.opticalLim);
        case 2 %Plot EnergyFlux Map
            pfisrData = get_2D_plot_inputs_at_time(inputH5FileStr,...
                'plotModeStr',expectedMaps{2},...
                'plotData',select_map_data(thisLayer,p.Results.map1Data,p.Results.map2Data,p.Results.map3Data),...
                'thisTimeIndx', thisTimeIndx(thisLayer));
           
            latWidth = p.Results.latLim(2)-p.Results.latLim(1);
            lonWidth = p.Results.lonLim(2)-p.Results.lonLim(1);
            
            pfisrData.diffEnergyFlux(imag(pfisrData.diffEnergyFlux(:))~=0)=nan;
            axesHandle(thisLayer)=axes;
            [axesmHandle(thisLayer),hEnergy]=plot_2D_energy_slice_geodetic_v2018...
                (pfisrData.diffEnergyFlux,...
            pfisrData.latitude, pfisrData.longitude,...
            pfisrData.zEnergyBin, pfisrData.thisTime,...
            p.Results.energySlice,latWidth,lonWidth,true,...
            p.Results.setEnergyTimeLabelOn);
            
            colormap(axesHandle(thisLayer),'inferno');
            cb(thisLayer) = colorbar(axesHandle(thisLayer),'westoutside');
            ylabel(cb(thisLayer),'log_1_0 [eV m^-^2 s^-^1 sr^-^1 eV^-^1]');
            caxis(axesHandle(thisLayer),p.Results.eFluxLim);
            
        case 3
            magFieldData = get_2D_plot_inputs_at_time(inputH5FileStr,...
                'plotModeStr',expectedMaps{3},...
                'plotData',select_map_data(thisLayer,p.Results.map1Data,p.Results.map2Data,p.Results.map3Data),...
                'thisTimeIndx', thisTimeIndx(thisLayer),...
                'magFieldModelStr',p.Results.magneticFieldModel);
            axesHandle(thisLayer)=axes;
            magLayerNo = thisLayer;
            
            if strcmp(p.Results.plotContours,'RE')
                 plotVariable = magFieldData.RE;
            elseif strcmp(p.Results.plotContours,'Lm')
                plotVariable = magFieldData.Lm;
            elseif strcmp(p.Results.plotContours,'Lstar')
                plotVariable = magFieldData.Lstar;
            end
            
            [axesmHandle(thisLayer), hMagnetic, cMagnetic] = plot_2D_magnetic_foot_points...
                (plotVariable, magFieldData.ionosphereCoord,...
                'BfieldModelStr',magFieldData.magFieldModelStr, 'thisTimeBfieldModel',magFieldData.thisTime,...
                'setMapOn',true, 'latLim', p.Results.latLim,'lonLim', p.Results.lonLim,'contourLineArray', p.Results.contourLineArray,...
                'setFieldLabelOn',false,'setTimeLabelOn', p.Results.setMagneticFieldTimeLabelOn);
        case 4
            warning(['No map in layer ',num2str(thisLayer)]);
    end
end

nAxes = length(axesHandle);
    if nAxes == 2
        linkaxes([axesHandle(1),axesHandle(2)]);
        set(axesHandle(2),'Position',get(axesHandle(1),'Position'),'XLim',get(axesHandle(1),'XLim'));
        if find(layer==3)>=1
            tempAxesHandle.XLim = axesHandle(magLayerNo).XLim;
        end
        setm(axesmHandle(2),'MapProjection','lambertstd',...
            'MapLatLimit',getm(axesmHandle(1),'MapLatLimit'),...
            'MapLonLimit',getm(axesmHandle(1),'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off');        
        setm(axesHandle(2),'FLatLimit',getm(axesmHandle(1),'FLatLimit'),...
            'FLonLimit',getm(axesmHandle(1),'FLonLimit'));
        alpha(axesHandle(2), 0.5);
        
    elseif nAxes == 3
        linkaxes([axesHandle(1),axesHandle(2),axesHandle(3)]);
        set(axesHandle(2),'Position',get(axesHandle(1),'Position'),'XLim',get(axesHandle(1),'XLim'));
        set(axesHandle(3),'Position',get(axesHandle(1),'Position'),'XLim',get(axesHandle(1),'XLim'));
        
        if find(layer==3)>=1
            tempAxesHandle.XLim = axesHandle(magLayerNo).XLim;
        end
        
        setm(axesmHandle(2),'MapProjection','lambertstd','MapLatLimit',getm(axesmHandle(1),'MapLatLimit'),...
            'MapLonLimit',getm(axesmHandle(1),'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off');

        setm(axesmHandle(3),'MapProjection','lambertstd','MapLatLimit',getm(axesmHandle(1),'MapLatLimit'),...
            'MapLonLimit',getm(axesmHandle(1),'MapLonLimit'),...
        'Frame','on','Grid','off','MeridianLabel','off','ParallelLabel','off',...
        'PLineVisible','off','MLineVisible','off');
        setm(axesHandle(2),'FLatLimit',getm(axesmHandle(1),'FLatLimit'),...
            'FLonLimit',getm(axesmHandle(1),'FLonLimit'));
        setm(axesHandle(3),'FLatLimit',getm(axesmHandle(1),'FLatLimit'),...
            'FLonLimit',getm(axesmHandle(1),'FLonLimit'));
                       
        alpha(axesHandle(2), 0.5);
        alpha(axesHandle(3), 0.5);

    end
    
        if layer(1)||layer(2)||layer(3)==3 && p.Results.setFieldLabelOn        
        htext = clabelm(cMagnetic,hMagnetic,p.Results.contourLabelArray,...
            'LabelSpacing',723);
        set(htext,'BackgroundColor','cyan','margin',2);
        end
        
        if find(layer==3)>=1
            set(axesHandle(magLayerNo),'XLim',tempAxesHandle.XLim);
        end
        

end

function OK = isfigure(h)
 if strcmp(get(h,'type'), 'figure')
     OK = true;
 else
     OK = false;
 end
     
end
 
function plotData = select_map_data(layer,map1Data,map2Data,map3Data)
    if layer==1
        plotData = map1Data;
    elseif layer==2
        plotData = map2Data;
    elseif layer==3
        plotData = map3Data;
    else
        plotData = struct();
    end
end