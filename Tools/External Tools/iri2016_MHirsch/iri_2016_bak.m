function data = iri_2016(time,altKm,glat,glon,setSamplePlot)
%IRI_2016 Matlab-Python wrapper of IRI2016 by MHirsch
%   time,glat,glon are one value each, while altKm can be a 1D array.

% https://www.scivision.co/matlab-python-user-module-import/
assert(~verLessThan('matlab', '9.5'), 'Matlab >= R2018b required, with python 3.6')

narginchk(4,5)
validateattributes(altKm, {'numeric'}, {'positive', 'vector'})
validateattributes(glat, {'numeric'}, {'scalar'})
validateattributes(glon, {'numeric'}, {'scalar'})

if nargin<5
    setSamplePlot = false;
end

validateattributes(setSamplePlot, {'logical'}, {'scalar'})

isoDateFormat = 'yyyy-mm-ddTHH:MM:SS';

switch class(time)
    case {'datetime', 'double'}, time = datestr(time, 30);
end

iono = py.iri2016.IRI(time, altKm, glat, glon);

data.Ne = xarray2mat(iono{'ne'});
data.Tn = xarray2mat(iono{'Tn'});
data.Ti = xarray2mat(iono{'Ti'});
data.Te = xarray2mat(iono{'Te'});

data.nOI = xarray2mat(iono{'nO+'});
data.nO2I = xarray2mat(iono{'nO2+'});
data.nHI = xarray2mat(iono{'nH+'});
data.nHeI = xarray2mat(iono{'nHe+'});
data.nNOI = xarray2mat(iono{'nNO+'});
data.nCI = xarray2mat(iono{'nCI'});
data.nNI = xarray2mat(iono{'nN+'});

data.totalIonDensity = (data.nOI + data.nCI +...
    data.nHeI + data.nHI + data.nNI + data.nNOI + data.nO2I);

data.NmF2 = xarray2mat(iono{'NmF2'});
data.hmF2 = xarray2mat(iono{'hmF2'});
data.NmF1 = xarray2mat(iono{'NmF1'});
data.hmF1 = xarray2mat(iono{'hmF1'});
data.NmE = xarray2mat(iono{'NmE'});
data.hmE = xarray2mat(iono{'hmE'});

data.B0 = xarray2mat(iono{'B0'});

data.coordinates.altKm = xarrayind2vector(iono, 'alt_km');
data.coordinates.lat = xarrayind2vector(iono, 'lat');
data.coordinates.lon = xarrayind2vector(iono, 'lon');

data.attribute.f107 = str2double(char(iono.attrs{'f107'}));
data.attribute.ap = str2double(char(iono.attrs{'ap'}));
data.attribute.time = datenum(char(iono.attrs{'time'}.isoformat()),isoDateFormat);

data.description = 'SI Units';

if setSamplePlot
    plot_iono(data);
end 

end

function plot_iono(data)

% p = create_panels(figure,'totalPanelNo',2,'panelHeight',60,...
%     'panelBreadth',60,'demargin',15,'marginleft',15,'marginright',15);
% q = p(1);
figure;
% q(1).select();
semilogx(data.Ne, data.coordinates.altKm)

hold on;
semilogx(data.nOI, data.coordinates.altKm)
hold on;
semilogx(data.nO2I, data.coordinates.altKm)
hold on;
semilogx(data.nNOI, data.coordinates.altKm)
hold on;
semilogx(data.nHI, data.coordinates.altKm)
hold on;
semilogx(data.nHeI, data.coordinates.altKm)
hold on;
semilogx(data.nCI, data.coordinates.altKm)
hold on;
semilogx(data.nNI, data.coordinates.altKm)
set(gca,'xscale','log')
legend('Ne','nO+','nO2+','nNO+','nH+','nHe+','nClusterIons','nN+','Location','northwest');
xlim([10^4 10^12]);

title({[datestr(data.attribute.time),' deg.  (',num2str(data.coordinates.lat),', ', num2str(data.coordinates.lon),')']})
xlabel('Density [m^-3]')
ylabel('altitude [km]')
grid('on')


end

% 
% function M = pyarray2mat(V)
% % M = double(py.array.array('d',py.numpy.nditer(py.numpy.asfortranarray(V)))); 
% M = double(py.numpy.asfortranarray(V));
% end

% function I = xarrayind2vector(V,key)
% 
% C = cell(V.indexes{key}.values.tolist);  % might be numeric or cell array of strings
% 
% if iscellstr(C) || isa(C{1}, 'py.str')
%     I = cellfun(@char, C, 'uniformoutput', false);
% elseif isa(C{1}, 'py.datetime.datetime')
%     I = char(C{1}.isoformat());
% else
%     I = cell2mat(C);
% end % if
% 
% end