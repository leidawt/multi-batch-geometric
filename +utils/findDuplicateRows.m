function [hasDuplicates, ixDupRows, dupRowValues] = findDuplicateRows(x)
    % find the duplicate rows of matrix x
    % hasDuplicates: whether there are duplicate rows
    % ixDupRows: index of duplicate rows
    % dupRowValues： duplicate rows
    % https://www.mathworks.com/matlabcentral/answers/184351-how-can-i-find-identical-rows-in-a-matrix
    [u, I, J] = unique(x, 'rows', 'first');
    hasDuplicates = size(u, 1) < size(x, 1);
    ixDupRows = setdiff(1:size(x, 1), I);
    dupRowValues = x(ixDupRows, :);
end
