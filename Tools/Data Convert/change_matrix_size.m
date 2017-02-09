function [ newMatrix ] = change_matrix_size( matrix, newRowSize, newColSize )
%% change_matrix_size.m Changes the size of the matrix by interpolating within
%--------------------------------------------------------------------------
% Input :
%-------
% matrix - Any real valued matrix
% newRowSize - Row size of the new matrix
% newColSize - Col size of the new matrix
%--------------------------------------------------------------------------
% Output :
%---------
% newMatrix - The new matrix interpolated with the new row and column size
%----------------------------------------------------------------------------
% Modified: 24th Jan 2017 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------

A = matrix;

oldRowSize = size(A,1);
oldColSize = size(A,2);

iOldRow = linspace(1,size(A,1),size(A,1));
iOldCol = linspace(1,size(A,2),size(A,2));

[X, Y] = meshgrid(iOldRow, iOldCol);

iNewRow = linspace(1,size(A,1),newRowSize);
iNewCol = linspace(1,size(A,2),newColSize);

[Xq, Yq] = meshgrid(iNewRow, iNewCol);

% Interpolation
newMatrix = interp2(X,Y,A,Xq,Yq,'nearest',0);



end

