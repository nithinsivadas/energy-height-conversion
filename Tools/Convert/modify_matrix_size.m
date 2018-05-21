function [ newMatrix ] = modify_matrix_size( matrix, newRowSize, newColSize )
%% modify_matrix_size.m Changes the size of the matrix by interpolating within
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
% Modified: 18th May 2018 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
% Notes   : 18th May 2018 - tried to improve the speed, and reduce
%           redundancy 
%----------------------------------------------------------------------------
if size(matrix,1) ~= newRowSize || size(matrix,2) ~= newColSize
    oldRowSize = size(matrix,1);
    oldColSize = size(matrix,2);

    iOldRow = linspace(1,oldRowSize,oldRowSize);
    iOldCol = linspace(1,oldColSize,oldColSize);

    [X, Y] = meshgrid(iOldRow, iOldCol);

    iNewRow = linspace(1,oldRowSize,newRowSize);
    iNewCol = linspace(1,oldColSize,newColSize);

    [Xq, Yq] = meshgrid(iNewRow, iNewCol);

    % Interpolation
    newMatrix = interp2(X,Y,matrix,Xq,Yq,'nearest',0);
else
    newMatrix = matrix;
end


end

