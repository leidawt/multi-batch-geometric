function C = cartesian(varargin)
    %cartesian product that can handle any type of input
    %https://stackoverflow.com/questions/9834254/cartesian-product-in-matlab
    %example cartesian([1 2],[3 4],[6 7])
    %ans =
    % 1     3     6
    % 1     3     7
    % 1     4     6
    % 1     4     7
    % 2     3     6
    % 2     3     7
    % 2     4     6
    % 2     4     7
    % or input a cell, where each cell element is a list

    n = length(varargin);

    if n == 1
        args = varargin{1};
        n = length(args);

        if length(args) == 1
            C = args{1};
            C = C';
            return;
        end

    else
        args = varargin;
    end

    [F{1:n}] = ndgrid(args{:});

    for i = n:-1:1
        G(:, i) = F{i}(:);
    end

    C = unique(G, 'rows');
end
