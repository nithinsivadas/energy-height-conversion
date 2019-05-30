function cmap = get_colormap3(lowerColor,midColor,upperColor,L)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if nargin<4
    L=32;
else 
    L = round(L/2);
end

if isstring(lowerColor)||ischar(lowerColor)
    rgb1 = get_rgb_triplet(lowerColor);
else
    rgb1 = lowerColor;
end

if isstring(midColor)||ischar(midColor)
    rgb2 = get_rgb_triplet(midColor);
else
    rgb2 = midColor;
end

if isstring(upperColor)||ischar(upperColor)
    rgb3 = get_rgb_triplet(upperColor);
else
    rgb3 = upperColor;
end

cmap0(:,1) = linspace(rgb1(1),rgb2(1),L);
cmap0(:,2) = linspace(rgb1(2),rgb2(2),L);
cmap0(:,3) = linspace(rgb1(3),rgb2(3),L);

cmap1(:,1) = linspace(rgb2(1),rgb3(1),L);
cmap1(:,2) = linspace(rgb2(2),rgb3(2),L);
cmap1(:,3) = linspace(rgb2(3),rgb3(3),L);

cmap = [cmap0;cmap1];
end

function rgb=get_rgb_triplet(colorStr)
    color = string(['y','m','c','r','g','b','w','k']');
    rgbArr = [1 1 0; 1 0 1;  0 1 1; 1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 0 0];    
    rgb = rgbArr(strcmp(color,string(colorStr)),:);
end