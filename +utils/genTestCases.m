function [arg_list, metaData] = genTestCases(Nlines, saveName, M, B_NUM, B_MIN, B_MAX, R_MIN, R_MAX, E_MIN, E_MAX, N_RANGE_MIN, N_RANGE_MAX)
    % generate test cases and save into to a .mat file
    % e.g. genTestCases(100,'testCase',5,10,7,15,0.2,0.4,0.7,0.95,2,4)
    if nargin < 3
        M = 5; % number of machines
    end

    if nargin < 4
        B_NUM = 10; % number of batches
    end

    if nargin < 5
        B_MIN = 7; % minimum batch size
    end

    if nargin < 6
        B_MAX = 15; % maximum batch size
    end

    if nargin < 7
        R_MIN = 0.2; % minimum repair rate
    end

    if nargin < 8
        R_MAX = 0.4; % maximum repair rate
    end

    if nargin < 9
        E_MIN = 0.7; % minimum machine efficiency
    end

    if nargin < 10
        E_MAX = 0.95; % maximum machine efficiency
    end

    if nargin < 11
        N_RANGE_MIN = 2; % minimum buffer size ratio
    end

    if nargin < 12
        N_RANGE_MAX = 4; % maximum buffer size ratio
    end

    metaData = struct( ...
        'N_LINES', Nlines, ...
        'M', M, ...
        'B_NUM', B_NUM, ...
        'B_MIN', B_MIN, ...
        'B_MAX', B_MAX, ...
        'R_MIN', R_MIN, ...
        'R_MAX', R_MAX, ...
        'E_MIN', E_MIN, ...
        'E_MAX', E_MAX, ...
        'N_RANGE_MIN', N_RANGE_MIN, ...
        'N_RANGE_MAX', N_RANGE_MAX ...
    );
    % record metadata
    arg_list = zeros(Nlines, 2 * M * B_NUM + M - 1 + B_NUM);

    for ind = 1:Nlines
        r = unifrnd(R_MIN, R_MAX, B_NUM, M);
        e = unifrnd(E_MIN, E_MAX, B_NUM, M);
        p = r ./ e - r;
        maxMTTR = max(1 ./ r);
        N = zeros(1, M - 1);

        for m = 1:M - 1
            N(m) = randi([N_RANGE_MIN, N_RANGE_MAX]) * fix(max(maxMTTR(m:m + 1)));
        end

        B = randi([B_MIN, B_MAX], 1, B_NUM);
        arg_list(ind, :) = utils.lineArgs2Vector(p, r, N, B);

    end

    if length(saveName) ~= 0
        % save
        save_str = sprintf("./%s_%d_%d.mat", saveName, M, B_NUM);
        save(save_str, 'arg_list', 'metaData');
    end

end
