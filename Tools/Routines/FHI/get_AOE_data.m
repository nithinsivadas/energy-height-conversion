function data=get_AOE_data(startTime, endTime)
    api = 'https://aoe2.net';
    startT = posixtime(startTime);
    endT   = posixtime(endTime);
    iTime = startT;
    k = 1;
    while iTime<endT
    url = [api '/api/matches' '?game=aoe2de&count=1000&since=',num2str(iTime)];
    [~,temp1] = system(['curl -sX GET "' url '"']);
    temp = jsondecode(temp1);
%   options = weboptions('ContentType','json');
%   temp = webread(url,options);
     l = length(temp);
     data(k:1:(k-1)+l) =  temp;
     k = k+l;
     iTime = data(end).finished;
     disp([num2str(100*(1-(endT-iTime)./(endT-startT)),'%10.2f') '%']);
    end
end