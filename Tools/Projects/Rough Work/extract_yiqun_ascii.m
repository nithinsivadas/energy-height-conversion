function [data] = extract_yiqun_ascii(fileStr)
%extract_yiqun_ascii Extract information from the specific ascii file
% mainly Hall & Pedersen conductivities from s/c loss-cone flux estimates. 
fileID = fopen(fileStr);
i = 1;
while ~feof(fileID)
    file(i) = string(fgetl(fileID));
    i = i+1;
end
i = 1;
nHeight=sscanf(file(1),' nHeight: %f'); 
nFile = length(file);
header1 = file(2:nHeight+2:nFile);
header2 = file(3:nHeight+2:nFile);
nTime = (nFile-1)./(nHeight+2);
timeIndx = 4:nHeight+2:nFile;
while i <= nTime
    j = 1;
    while j <= nHeight
        temp = sscanf(file(timeIndx(i)+j-1),'%f');
        arr(i,j,:)=temp(2:end);
        if i==1
            height(j) = temp(1);
        end
        j = j+1;
    end
    data.header1(i,:) = sscanf(header1(i),'%f');
    i = i+1;
end
data.ionRate = squeeze(arr(:,:,1));
data.NeIn = squeeze(arr(:,:,2));
data.NeOut = squeeze(arr(:,:,3));
data.nOI = squeeze(arr(:,:,4));
data.nO2I = squeeze(arr(:,:,5));
data.nNOI = squeeze(arr(:,:,6));
data.sigmaP = squeeze(arr(:,:,7));
data.sigmaH = squeeze(arr(:,:,8));
data.alt = height';
data.time = datenum(data.header1(:,1:6));
data.lat = data.header1(:,7);
data.lon = data.header1(:,8);
data.mlat = data.header1(:,9);
data.mlon = data.header1(:,10);
data.mlt = data.header1(:,11);
fclose(fileID);

end

