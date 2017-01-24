function [ rotatedArray ] = rotate_array( array, alpha )
%Rotate the array clockwise w.r.t data

a1 = array - alpha;
I = find(a1<0);
a1(I) = 360 + a1(I);
rotatedArray = a1;
end

