function el1 = tilt_elevation(el,del)
%tile_elevation This function can rotate elevation angle by preserving the
%wrapping
%   Detailed explanation goes here
el1 = 90.*sawtooth(deg2rad(el + del + 90),1/2);

end

