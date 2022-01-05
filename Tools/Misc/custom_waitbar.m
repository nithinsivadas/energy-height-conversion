function custom_waitbar(hWaitBar,i,n,commentStr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
waitbar(i/n,hWaitBar,[commentStr,' ',sprintf('%3.0f',100*i/n),'%'])
end

