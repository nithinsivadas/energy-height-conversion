function cmap = get_colormap(lowerColor,upperColor)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
rgb1 = get_rgb_triplet(lowerColor);
rgb2 = get_rgb_triplet(upperColor);
L=64;
cmap(:,1) = linspace(rgb1(1),rgb2(1),L);
cmap(:,2) = linspace(rgb1(2),rgb2(2),L);
cmap(:,3) = linspace(rgb1(3),rgb2(3),L);

end

function rgb=get_rgb_triplet(colorStr)
    color = string(['y','m','c','r','g','b','w','k']');
    rgbArr = [1 1 0; 1 0 1;  0 1 1; 1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 0 0];    
    rgb = rgbArr(strcmp(color,string(colorStr)),:);
end