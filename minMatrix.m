function [row,col] = minMatrix(A)
    minMatrix = min(A(:));
    [row,col] = find(A==minMatrix);
end