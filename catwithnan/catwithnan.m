function new = catwithnan(cellArray, dim)
%CATWITHNAN Concatenate arrays into NaN-delimited array
%
% new = catwithnan(cellArray, dim)
%
% Concatenates a cell array into a matrix along the specified dimension,
% separating former cells with NaNs.
%
% Input variables:
%
%   cellArray:  cell array to be concatenated
%
%   dim:        1 or 2, dimension along which cellArray will be
%               concatenated.  If dim = 1, all elements in the cell array
%               must have the same number of columns; if dim = 2, all
%               elements must have the same number of rows.
%
% Output variables:
%
%   new:        concatenated array
%
% Example:
%
% a = catwithnan({ones(2,3), rand(5,3), [1 2 3]}, 1)
% 
% a =
% 
%     1.0000    1.0000    1.0000
%     1.0000    1.0000    1.0000
%        NaN       NaN       NaN
%     0.4057    0.0579    0.2028
%     0.9355    0.3529    0.1987
%     0.9169    0.8132    0.6038
%     0.4103    0.0099    0.2722
%     0.8936    0.1389    0.1988
%        NaN       NaN       NaN
%     1.0000    2.0000    3.0000
%        NaN       NaN       NaN

% Copyright 2005 Kelly Kearney

ncells = numel(cellArray);
new = [];

if dim == 1
    cellArray = cellArray(:);
    nrows = cellfun('size', cellArray, 1);
    ncols = cellfun('size', cellArray, 2);
    ncols = unique(ncols);
    if ~isscalar(ncols)
        error('All arrays in the cell array must have the same number of columns');
    end
    startIndex = 1;
    new = NaN(sum(nrows) + numel(cellArray), ncols);
    
    for icell = 1:length(cellArray)
        endIndex = startIndex + nrows(icell) - 1;
        new(startIndex:endIndex, :) = cellArray{icell};
        startIndex = endIndex + 2;
    end
elseif dim == 2
    cellArray = cellArray(:)';
    nrows = cellfun('size', cellArray, 1);
    ncols = cellfun('size', cellArray, 2);
    nrows = unique(nrows);
    if ~isscalar(nrows)
        error('All arrays in the cell array must have the same number of rows');
    end
    startIndex = 1;
    new = NaN(nrows, sum(ncols) + numel(cellArray));
    for icell = 1:length(cellArray)
        endIndex = startIndex + ncols(icell) - 1;
        new(:, startIndex:endIndex) = cellArray{icell};
        startIndex = endIndex + 2;
    end
end
        
        
