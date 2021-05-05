%% Error estimation from solar wind data

dataFolder = 'G:\My Drive\Research\Projects\Data\omni_components\';

ace1 = extract_set_of_data(dataFolder,'ace');
wind1 = extract_set_of_data(dataFolder,'wind');
% dscov = extract_data('G:\My Drive\Research\Projects\Data\omni_components\dscov_min_b2017.txt');
% wind = extract_data('G:\My Drive\Research\Projects\Data\omni_components\wind_min_b2017.txt');

%%
[~,ai,wi]=intersect(ace1.datetime,wind1.datetime);
ace = ace1(ai,:);
wind = wind1(wi,:);

%%

%% Specifying a grid
gride = -5:0.1:5;
gridx = 0:0.4:40;
[E, X] = meshgrid(gride,gridx);

% Generating a Kernal density estimate, to have a model of the continous
% probability distribution of error

[fx,xi] = ksdensity((wind.E_kl*10^-3),gridx);
Fx = griddedInterpolant(xi,fx);
[fex,xii] = ksdensity([(ace.E_kl-wind.E_kl)*10^-3, ace.E_kl*10^-3],[E(:),X(:)]);
Fex = scatteredInterpolant(xii(:,1),xii(:,2),fex);

Dx = diff(gridx);
Dx(end+1) = Dx(end);

% [fx,xi] = ksdensity((wind.E_kl*10^-3));
% Fx = griddedInterpolant(xi,fx);
% [fex,xii] = ksdensity([(ace.E_kl-wind.E_kl)*10^-3, ace.E_kl*10^-3]);
% Fex = scatteredInterpolant(xii(:,1),xii(:,2),fex);


gride1 = -5:0.1:5;
gridx1 = 0:0.1:40;
[E1, X1] = meshgrid(gride1,gridx1);
figure; 
p=pcolor(E1,X1,(Fex(E1,X1)./(Fx(gridx1)')));
set(p, 'EdgeColor', 'none');
colorbar;
caxis([0,2]);
xlim([-5,5]);
colormap('inferno');
xlabel('p( \epsilon ,E_{M} ) [mV/m]');
ylabel('E_{M} [mV/m]');

%%
% T = extract_set_of_data(dataFolder,satellite);
% wind = extract_set_of_data(dataFolder,'wind');

function T = extract_set_of_data(loadFolder,satellite)
    fileStr = get_files_in_folder(loadFolder,[satellite,'_min_b*.txt']);
    for i=1:1:length(fileStr)
        T1 = extract_data([loadFolder,fileStr{i}]);
        if i==1
            T=T1;
        else
            T = [T;T1];
        end
    end
end


function T1 = extract_data(loadFile)
 
        format = '%4f %4f %3f %3f %4f %4f %4.1f %6f %6.2f %6.2f %6.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7f %6.2f %8.2f %8.2f %4f %8.1f %8.1f %8.1f %8.1f %7.2f %9.0f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7f %7f';
        tempData = load_ascii_files(loadFile, format, 0);
        T.datetime = datetime(tempData{1},1,tempData{2},tempData{3},tempData{4},zeros(size(tempData{4})));
        T.timeshift = tempData{8}; T.timeshift(T.timeshift==999999) = nan;
        T.BxGSM = tempData{13};  T.BxGSM(T.BxGSM > 9999) = nan;
        T.ByGSM = tempData{16}; T.ByGSM(T.ByGSM > 9999) = nan;
        T.BzGSM = tempData{17}; T.BzGSM(T.BzGSM > 9999) = nan;
        T.velocity = tempData{23}; T.velocity(T.velocity > 99999) = nan;
        T.B_T = (T.ByGSM.^2 + T.BzGSM.^2).^0.5;
        T.theta_kl = wrapTo2Pi(atan2(T.B_T,T.BzGSM));
        T.E_kl = T.velocity.*T.B_T.*(sin(T.theta_kl/2)).^2;
        T.XGSE = tempData{29}; T.XGSE(T.XGSE > 9999) = nan;
        T.YGSE = tempData{30}; T.YGSE(T.YGSE > 9999) = nan;
        T.ZGSE = tempData{31}; T.ZGSE(T.ZGSE > 9999) = nan;
        T.noseXGSE = tempData{32}; T.noseXGSE(T.noseXGSE > 9999) = nan;
        T.noseYGSE = tempData{33}; T.noseYGSE(T.noseYGSE > 9999) = nan;
        T.noseZGSE = tempData{34}; T.noseZGSE(T.noseZGSE > 9999) = nan;
        T.distance = sqrt((T.XGSE-T.noseXGSE).^2 + (T.YGSE-T.noseYGSE).^2 + (T.ZGSE-T.noseZGSE).^2);
        T1 = struct2table(T);
end
