clear variables;
iTime = 1:20;
% integrationTimeArray =linspace(500,5000,length(iTime));
nPulseTrainsArray = iTime;

iSNR = 1:20;
SNRArray = logspace(-1,log10(2),length(iSNR));

ifWidth = 1:10;
fWidthArray = linspace(1,100,length(ifWidth));

mean_dist=@(dist) mean(real(dist(:)));
std_dist=@(dist) std(real(dist(:)));
multiWaitbar('CLOSEALL');
multiWaitbar('Integration time',0);
multiWaitbar('SNR',0);

% multiWaitbar('Spectral width',0);
tic
for t=iTime
    multiWaitbar('Integration time','increment',1/iTime(end),'color','g');
    for s=iSNR
    multiWaitbar('SNR','increment',1/iSNR(end));
        for fw=ifWidth

        % Controllable Inputs
        nExperiments = 2^8;
%         integrationTime = integrationTimeArray(t); % ms
        
        % ACF Characteristics
        SNR = SNRArray(s);
%         altitude = 60; %km 
        
        noisePower = 100;
        signalPower = SNR*noisePower;
        fWidth = fWidthArray(fw); % in Hz
        dopplerFrequency = 10; %in Hz
        
        
        
        % Pulse characteristics
        IPP = 2; %in ms (inter pulse period)
%         nyqyustSamplingRate = 1/(IPP*10^-3);
%         if IPP ==0
%             nyquistSamplingRate = 10*2*fWidth; % n times the nyquist frequency
%         elseif (IPP*10^-3<(1/2*fWidth)) % if IPP is less than nyquist sampling time
%             nyquistSamplingRate=1/(IPP*10^-3);
%         else
%             warning('IPP is longer than 1/(2*fWidth) - the signal is being undersampled');
%         end
        
        maxLag = 4*(1/fWidth); % n times the decorrelation time scale
        
        % Sampling
%         sampleTime = 1/nyquistSamplingRate; %IPP
        sampleTime = IPP*10^-3;
        nSamples = 2^nextpow2(maxLag/sampleTime);
%         nPulseTrains = ceil(integrationTime/(nSamples*sampleTime*10^3));
        nPulseTrains = nPulseTrainsArray(t);
        
        
%         display(['Sample Time/ IPP: ',num2str(sampleTime*10^3,4),' ms']);
%         display(['Decorrelation Time: ',num2str((1/fWidth).*10^3,4),' ms']);
%         display(['Max Lag: ',num2str(nSamples*sampleTime*10^3,4),' ms']);
%         display(['nSamples: ',num2str(nSamples)]);
%         display(['Experiments #: ',num2str(nExperiments),' Pulse trains #: ',num2str(nPulseTrains)]);
%         display(['Integration time: ',num2str(nSamples*nPulseTrains*sampleTime*10^3,5),' ms']);
         
         [noiseDist,signalDist,dopplerShiftDist,fWidthDist]...
            = get_ACF_parameter_distribution_1(nSamples,nPulseTrains,nExperiments,...
            noisePower,signalPower,dopplerFrequency,fWidth,sampleTime);
         
         data.mean.noise(t,s,fw)=mean_dist(noiseDist);
         data.mean.signal(t,s,fw)=mean_dist(signalDist);
         data.mean.doppler(t,s,fw)=mean_dist(dopplerShiftDist);
         data.mean.fWidth(t,s,fw)=mean_dist(fWidthDist);
         data.std.noise(t,s,fw)=std_dist(noiseDist);
         data.std.signal(t,s,fw)=std_dist(signalDist);
         data.std.doppler(t,s,fw)=std_dist(dopplerShiftDist);
         data.std.fWidth(t,s,fw)=std_dist(fWidthDist);
         data.exp.IPP(t,s,fw)=sampleTime;
         data.exp.decorrelationTime(t,s,fw)=1/fWidth;
         data.exp.maxLag(t,s,fw) = nSamples*sampleTime;
         data.exp.SNR(t,s,fw) = SNR;
         data.exp.integrationTime(t,s,fw) =nSamples*nPulseTrains*sampleTime;
         data.exp.nExperiments(t,s,fw) =nExperiments;
         data.exp.nPulseTrains(t,s,fw) =nPulseTrains;

        end
%         multiWaitbar('Spectral width','reset');
    end
    multiWaitbar('SNR','reset');
end

toc

% Checking variations of different parameter error
data.mean.SNR = data.mean.signal./data.mean.noise;
data.fracError.signal = data.std.signal./data.mean.signal;
data.fracError.noise = data.std.noise./data.mean.noise;
data.fracError.SNR = data.fracError.signal+abs(data.fracError.noise);

% Store the variable
storeDir = [initialize_root_path,'LargeFiles',filesep,'ISRMonteCarlo',filesep];
fileNo=1;
storeFileName=['data',num2str(fileNo),'.mat'];
while isfile([storeDir,storeFileName])
    fileNo=fileNo+1;
    storeFileName=['data',num2str(fileNo),'.mat'];
end
save([storeDir,storeFileName],'data');


    
