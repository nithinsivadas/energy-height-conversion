function maginput = filter_irbem_maginput(magFieldNo,maginput)
%filter_irbem_maginput Filter irbem magnetic field inputs to appropriate limits.
% Only works for magfieldNo 6,7(T96),9(T01)
%   Detailed explanation goes here

if magFieldNo == 6 % Olson & Pfitzer dynamic
   maginput(:,3) = out_of_range(maginput(:,3),5,50); %density
   maginput(:,4) = out_of_range(maginput(:,4),300,500); %velocity
   maginput(:,2) = out_of_range(maginput(:,2),-100,20); % Dst

elseif magFieldNo == 7 %TS1996
   maginput(:,2) = out_of_range(maginput(:,2),-100,20); % Dst
   maginput(:,5) = out_of_range(maginput(:,5),0.5,10); % Pressure
   maginput(:,6) = out_of_range(maginput(:,6),-10,10); % ByIMF
   maginput(:,7) = out_of_range(maginput(:,7),-10,10); % BzIMF
elseif magFieldNo == 9 %TS2001
   maginput(:,2) = out_of_range(maginput(:,2),-50,20); % Dst
   maginput(:,5) = out_of_range(maginput(:,5),0.5,5); % Pressure
   maginput(:,6) = out_of_range(maginput(:,6),-5,5); % ByIMF
   maginput(:,7) = out_of_range(maginput(:,7),-5,5); % BzIMF
   maginput(:,8) = out_of_range(maginput(:,8),0,10); % G1
   maginput(:,9) = out_of_range(maginput(:,9),0,10); % G2
else
  warning('This version of filter_irbem_magnipnut does not filter inputs for this magnetic field model');
end


end

function x = out_of_range(x,xmin,xmax)
    x(x<xmin)=xmin;
    x(x>xmax)= xmax;
end
