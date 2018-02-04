function [X,Y,Z,P] = reshapegen(coords,p)
% reshapegen.m
% by John Swoboda
% This function is used to reshape data that may have not been flattened 
% using MATLAB's ordering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% coords - A nx3 array locations where it is assumed that thos locations
% come from flattened three-d arrays.
% p - A nx1 array of data points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% X Y Z - Two or Three dimensional arrays formed from the mesh grid of unique
%  values in each column in coords.
% P - The two or three dimensional array that is the reshaped data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xvec,~,icx] = unique(coords(:,1));
[yvec,~,icy] = unique(coords(:,2));
[zvec,~,icz] = unique(coords(:,3));

[X,Y,Z] = meshgrid(xvec,yvec,zvec);
P = zeros(size(X));
linearInd = sub2ind(size(X), icy,icx,icz);
P(linearInd) = p(:);



