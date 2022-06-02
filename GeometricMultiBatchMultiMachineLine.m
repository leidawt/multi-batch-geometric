classdef GeometricMultiBatchMultiMachineLine
    % GeometricMultiBatchMultiMachineLine

    properties
        M % number of machines
        B % batch size (vector of size 1*K)
        p % machine repair rate (matrix of size K*M)
        r % machine failure rate (matrix of size K*M)
        n % buffer size (vector of size 1*M-1)
        batchNum % number of batches
        totalNum % number of total parts
        efficienceMappingTable1 % lookup table for p
        efficienceMappingTable2 % lookup tabel for r
    end

    methods

        function obj = GeometricMultiBatchMultiMachineLine(p, r, n, B)
            % constructor
            %   p: machine failure rate (matrix of size K*M)
            %   r: machine repair rate (matrix of size K*M)
            %   n: buffer size (vector of size 1*M-1)
            %   B: batch size (vector of size 1*K)
            obj.M = size(p, 2);
            obj.batchNum = length(B);
            obj.totalNum = sum(B);
            obj.p = p;
            obj.n = n;
            obj.r = r;
            obj.B = B;
            assert(obj.batchNum == size(p, 1));
            assert(size(p, 2) == size(r, 2));
            assert(size(p, 1) == size(r, 1));
            assert(size(p, 1) == obj.batchNum);
            assert(length(n) + 1 == obj.M);
            % prepare lookup table for p, r
            obj.efficienceMappingTable1 = repelem(obj.p, obj.B, 1);
            obj.efficienceMappingTable1(end + 1, :) = obj.efficienceMappingTable1(end, :);
            obj.efficienceMappingTable1 = obj.efficienceMappingTable1';
            obj.efficienceMappingTable2 = repelem(obj.r, obj.B, 1);
            obj.efficienceMappingTable2(end + 1, :) = obj.efficienceMappingTable2(end, :);
            obj.efficienceMappingTable2 = obj.efficienceMappingTable2';

        end

        function [A, C] = calTransMatrix(obj, hFirst)
            % Calculate the one-step transition matrix A and the observation matrix C
            % C = [PR; CR; Pct]
            buffer_cumprod = cumprod([obj.n + 1, 1], 'reverse'); %if n=[2 2 3],then buffer_cumprod=[36 12 4 1]
            mList = {}; mList{obj.M - 1} = [];

            for m = 1:(obj.M - 1)
                mList{m} = 0:obj.n(m);
            end

            buffer_list = utils.cartesian(mList); % All possible buffer occupancy combinations

            args = struct();
            args.buffer_cumprod = buffer_cumprod(2:end);
            args.sMachine = 2^obj.M; % Number of machine state combinations
            args.sBuffer = buffer_cumprod(1); % Number of buffer state combinations
            args.M = obj.M;
            args.buffer_list = buffer_list;
            args.S = 2^obj.M * buffer_cumprod(1) * (1 + obj.totalNum);
            args.n = obj.n;
            args.efficienceMappingTable1 = obj.efficienceMappingTable1;
            args.efficienceMappingTable2 = obj.efficienceMappingTable2;
            args.totalNum = obj.totalNum;

            % for idx = 1:args.S
            %     state = idx2state(args, idx);
            %     idx_ = state2idx(args, state);
            %     disp([idx, idx_, state]);
            %     assert(idx == idx_);
            % end
            % stateList = zeros(args.S, (args.M * 2));

            % for idx = 1:args.S
            %     state = idx2state(args, idx);
            %     stateList(idx, :) = state;
            % end

            function idx = state2idx(args, state)
                % state=(h1,h2,...,s1,s2...f) -> index
                s_ = state(args.M:end - 1);
                h_ = state(1:args.M - 1);
                f_ = state(end);
                idx = 1 + s_ * (cumprod([2 * ones(1, args.M - 1) 1], 'reverse'))';
                idx = idx + (h_ * args.buffer_cumprod') * args.sMachine;
                idx = idx + f_ * (args.sMachine * args.sBuffer);
            end

            function state = idx2state(args, idx)
                % index -> state=(h1,h2,...,s1,s2...f)
                f_ = fix(idx / (args.sMachine * args.sBuffer));
                idx = rem(idx, (args.sMachine * args.sBuffer));

                if idx == 0
                    idx = args.sMachine * args.sBuffer;
                    f_ = f_ - 1;
                end

                mc_idx = rem(idx, args.sMachine);

                if mc_idx == 0
                    mc_idx = args.sMachine - 1;
                else
                    mc_idx = mc_idx - 1;
                end

                s_ = dec2bin(mc_idx, args.M) - '0';
                buffer_idx = floor((idx - 1) / args.sMachine);
                h_ = args.buffer_list(buffer_idx + 1, :);
                state = [h_, s_, f_];
            end

            function prob = cal_prob(args, s_old, s_new, h, f)
                % Calculate transition probability
                prob = 1;
                fs = [cumsum(h, 'reverse') + f, f];

                for m_ = 1:args.M
                    p_ = args.efficienceMappingTable1(m_, fs(m_) + 1);
                    r_ = args.efficienceMappingTable2(m_, fs(m_) + 1);

                    prob_temp = (s_old(m_) == 1 && s_new(m_) == 0) * p_ + ...
                        (s_old(m_) == 1 && s_new(m_) == 1) * (1 - p_) + ...
                        (s_old(m_) == 0 && s_new(m_) == 1) * r_ + ...
                        (s_old(m_) == 0 && s_new(m_) == 0) * (1 - r_);
                    prob = prob * prob_temp;

                end

            end

            function h_new = cal_h(args, h, f, ms)
                % update new h
                fs = [cumsum(h, 'reverse') + f, f];
                ms = ms .* (fs < args.totalNum);
                h_ = zeros(1, args.M - 1);
                h_new = zeros(1, args.M - 1);
                h_(end) = h(end) - ms(args.M) * min([h(end), 1]);

                for ii = args.M - 2:-1:1
                    h_(ii) = h(ii) - ms(ii + 1) * min([h(ii), args.n(ii + 1) - h_(ii + 1), 1]);
                end

                for ii = 2:args.M - 1
                    h_new(ii) = h_(ii) + ms(ii) * min([h(ii - 1), args.n(ii) - h_(ii), 1]);
                end

                h_new(1) = h_(1) + ms(1) * min([args.n(1) - h_(1), 1]);
            end

            A = sparse(args.S, args.S);
            C = sparse(3, args.S);

            for j = 1:args.S
                state = idx2state(args, j);
                f = state(end);
                s = state(args.M:end - 1);
                h = state(1:args.M - 1);
                isValidState = (sum(h) + f) <= args.totalNum;
                isProd = h(end) > 0 && s(end) && f < args.totalNum;

                if isValidState
                    f_new = f + isProd;
                    mStates = dec2bin(2^args.M - 1:-1:0) - '0';

                    for imStates = 1:length(mStates)
                        s_new = mStates(imStates, :);

                        if hFirst
                            h_new = cal_h(args, h, f, s); % update h first
                        else
                            h_new = cal_h(args, h, f, s_new); % update s first
                        end

                        i = state2idx(args, [h_new, s_new, f_new]);
                        A(i, j) = cal_prob(args, s, s_new, h, f);
                    end

                end

                C(1, j) = isProd;
                C(2, j) = (f + sum(h)) < args.totalNum && s(1) && (h(1) ~= args.n(1)); % this has a negligible error, should not be a problem
                C(3, j) = isProd && (f == (args.totalNum - 1));

            end

        end

        function [PR, CR, CT] = markovAnalysis(obj, A, C, SIM_T)
            % Exact analysis with Markov analysis
            PR = zeros(1, SIM_T);
            CR = zeros(1, SIM_T);
            CT = 0;
            x = sparse(size(A, 1), 1);
            x(1) = 1;

            for t = 1:SIM_T
                PR(t) = C(1, :) * x;
                CR(t) = C(2, :) * x;
                CT = CT + t * C(3, :) * x;
                x = A * x;
            end

        end

        function [PR, CR, WIP, ST, BL, CT] = MonteCarloAnalysis(obj, SIM_T, ITER, h_0, s_0)
            % Monte Carlo analysis of geometric lines
            % SIM_T: number of evaluation steps
            % ITER: number of repetitions
            % h_0: initial buffer occupancy
            % s_0: initial machine status

            function [p, r] = getMachineEfficience(m, count)
                %getMachineEfficience - Returns the efficiency of the part based on the entered material count value
                %
                % Syntax: [p, r] = getMachineEfficience(m, count)
                % Accelerate this function with cache
                if count <= obj.totalNum
                    p = obj.efficienceMappingTable1(m, count);
                    r = obj.efficienceMappingTable2(m, count);
                else
                    p = 0;
                    r = 0;
                end

            end

            % Build virtual production line to facilitate simulation calculation
            % count is the number of parts that have been processed by the machine
            line = repmat(struct('p', 0, 'r', 0, 'n', 0, 'buffer_occ', 0, ...
            'isBlock', false, 'isStarve', false, 'mState', false, 'count', 0), obj.M, 1);

            for m = 1:obj.M - 1
                line(m).p = obj.p(m);
                line(m).r = obj.r(m);
                line(m).n = obj.n(m);

                if nargin > 4
                    % use the given initial value
                    line(m).buffer_occ = h_0(m);
                    line(m).mState = s_0(m) == 1;
                end

            end

            % m=M
            line(obj.M).p = obj.p(obj.M);
            line(obj.M).r = obj.r(obj.M);

            if nargin > 4
                % use the given initial value
                line(obj.M).mState = s_0(obj.M) == 1;
            end

            % pre-allocate the result vector
            PR = zeros(1, SIM_T);
            CR = zeros(1, SIM_T);
            WIP = zeros(obj.M - 1, SIM_T);
            ST = zeros(obj.M, SIM_T);
            BL = zeros(obj.M, SIM_T);
            CT = zeros(obj.batchNum, 1);

            for iITER = 1:ITER
                line_ = line;
                complete_flag = zeros(obj.batchNum, 1); % batch completion flag

                for t = 1:SIM_T
                    p_rand = rand(obj.M, 1); % sampling
                    r_rand = rand(obj.M, 1); % sampling

                    for m = 1:obj.M
                        [p_, r_] = getMachineEfficience(m, line_(m).count + 1);
                        % determine machine switch status
                        line_(m).mState = p_rand(m) > p_ && line_(m).mState ...
                        || r_rand(m) < r_ && ~line_(m).mState;
                        % determine machine starvation status
                        if m == 1
                            line_(m).isStarve = false; % first machine never starve
                        else
                            line_(m).isStarve = line_(m - 1).buffer_occ == 0;
                        end

                    end

                    % determine machine blockage status
                    for m = obj.M - 1:-1:1
                        line_(m).isBlock = (line_(m).buffer_occ == line_(m).n) ...
                            && (~(line_(m + 1).mState && ~(line_(m + 1).isBlock)));
                    end

                    line_(obj.M).isBlock = false; % final machine never block

                    for m = 1:obj.M
                        % record BL ST
                        if m == 1
                            BL(m, t) = BL(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isBlock);
                        elseif m == obj.M
                            ST(m, t) = ST(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isStarve);
                        else
                            BL(m, t) = BL(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isBlock);
                            ST(m, t) = ST(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isStarve);
                        end

                        % determine production status, update WIP
                        isProd = (line_(m).mState) ...
                        && (line_(m).count < obj.totalNum) ...
                            && ~(line_(m).isBlock) && ~(line_(m).isStarve);

                        if isProd
                            line_(m).count = line_(m).count + 1;

                            if m == 1
                                line_(m).buffer_occ = line_(m).buffer_occ + 1;
                                CR(t) = CR(t) + 1;
                            elseif m == obj.M
                                line_(m - 1).buffer_occ = line_(m - 1).buffer_occ - 1;
                                PR(t) = PR(t) + 1;
                            else
                                line_(m).buffer_occ = line_(m).buffer_occ + 1;
                                line_(m - 1).buffer_occ = line_(m - 1).buffer_occ - 1;
                            end

                        end

                    end

                    % record WIP
                    for m = 1:obj.M - 1
                        WIP(m, t) = WIP(m, t) + line_(m).buffer_occ;
                    end

                    % record CT
                    c = cumsum(obj.B);

                    for k = 1:obj.batchNum

                        if line_(obj.M).count == c(k) && complete_flag(k) == 0
                            CT(k) = CT(k) + t;
                            complete_flag(k) = 1;
                        end

                    end

                end

            end

            % calculate sample expectation
            PR = PR ./ ITER;
            CR = CR ./ ITER;
            WIP = WIP ./ ITER;
            ST = ST ./ ITER;
            BL = BL ./ ITER;
            CT = CT ./ ITER;
        end

        function [PR, CR, WIP, ST, BL, CT] = MonteCarloAnalysisNonGeo(obj, SIM_T, ITER, h_0, s_0, beta)
            % Monte Carlo analysis for DW-distribution lines
            % SIM_T: number of evaluation steps
            % ITER: number of repetitions
            % h_0: initial buffer occupancy
            % s_0: initial machine status
            % beta shape parameter of the discrete Weibull distribution

            % prepare lookup table
            DWpmf_meoized = memoize(@utils.DWpmf);
            DWlookTable = zeros(1000, 2); %[p mean]
            ind = 1;

            for pDW = 0.001:0.001:1

                [~, mean, ~, ~] = DWpmf_meoized(pDW, beta);
                DWlookTable(ind, :) = [pDW, mean];
                ind = ind + 1;

            end

            function pmf = getDWpmf(mean)
                % get the pmf of the DW distribution given the mean
                assert(mean <= max(DWlookTable(:, 2)));
                assert(mean >= min(DWlookTable(:, 2)));
                [~, idx] = min(abs(DWlookTable(:, 2) - mean));
                pp = DWlookTable(idx, 1);
                pmf = DWpmf_meoized(pp, beta);
            end

            function [p, r] = getMachineEfficience(m, count)
                p = obj.efficienceMappingTable1(m, count);
                r = obj.efficienceMappingTable2(m, count);

            end

            batchED = cumsum(obj.B);
            batchST = batchED - obj.B + 1;

            function BatchType = getBatchType(count)
                %  determine the type of the workpiece given the production count

                flag = batchST <= count & batchED >= count;

                if sum(flag) == 0
                    % must be the final part
                    BatchType = length(batchED);
                else
                    BatchType = (1:obj.batchNum) * flag';
                end

            end

            line = repmat(struct('type', 1, 'n', 0, 'buffer_occ', 0, ...
                'isBlock', false, 'isStarve', false, 'mState', false, 'count', 0, 'T', 0), obj.M, 1);

            for m = 1:obj.M

                if nargin > 4
                    % use the given initial value
                    line(m).mState = s_0(m) == 1;
                end

                if m ~= obj.M

                    if nargin > 4
                        % use the given initial value
                        line(m).buffer_occ = h_0(m);
                    end

                    line(m).n = obj.n(m);
                end

            end

            % pre-allocate the result vector
            PR = zeros(1, SIM_T);
            CR = zeros(1, SIM_T);
            WIP = zeros(obj.M - 1, SIM_T);
            ST = zeros(obj.M, SIM_T);
            BL = zeros(obj.M, SIM_T);
            CT = zeros(obj.batchNum, 1);

            for iITER = 1:ITER
                line_ = line;

                for m = 1:obj.M
                    [p_, r_] = getMachineEfficience(m, line_(m).count + 1);
                    p_pmf = getDWpmf(1 / p_);
                    r_pmf = getDWpmf(1 / r_);

                    if line_(m).mState
                        % p
                        line_(m).T = randsample(1:length(p_pmf), 1, true, p_pmf);
                    else
                        % r
                        line_(m).T = randsample(1:length(r_pmf), 1, true, r_pmf);
                    end

                end

                complete_flag = zeros(obj.batchNum, 1); % batch completion flag

                for t = 1:SIM_T

                    for m = 1:obj.M
                        % determine machine work status
                        line_(m).T = line_(m).T - 1;

                        if line_(m).T == 0 || (getBatchType(line_(m).count + 1) ~= line_(m).type)

                            line_(m).type = getBatchType(line_(m).count + 1);

                            if line_(m).T == 0
                                line_(m).mState = ~line_(m).mState;
                            end

                            % sample a new T
                            [p_, r_] = getMachineEfficience(m, line_(m).count + 1);
                            p_pmf = getDWpmf(1 / p_);
                            r_pmf = getDWpmf(1 / r_);

                            if line_(m).mState
                                % p
                                line_(m).T = randsample(1:length(p_pmf), 1, true, p_pmf);
                            else
                                % r
                                line_(m).T = randsample(1:length(r_pmf), 1, true, r_pmf);
                            end

                        end

                        % determine machine starvation status
                        if m == 1
                            line_(m).isStarve = false; % first machine never starve
                        else
                            line_(m).isStarve = line_(m - 1).buffer_occ == 0;
                        end

                    end

                    % determine machine blockage status
                    for m = obj.M - 1:-1:1
                        line_(m).isBlock = (line_(m).buffer_occ == line_(m).n) ...
                            && (~(line_(m + 1).mState && ~(line_(m + 1).isBlock)));
                    end

                    line_(obj.M).isBlock = false; % final machine never block

                    for m = 1:obj.M
                        % record BL ST
                        if m == 1
                            BL(m, t) = BL(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isBlock);
                        elseif m == obj.M
                            ST(m, t) = ST(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isStarve);
                        else
                            BL(m, t) = BL(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isBlock);
                            ST(m, t) = ST(m, t) + 1 * (line_(m).count < obj.totalNum) * (line_(m).mState) * (line_(m).isStarve);
                        end

                        % determine production status and update WIP
                        isProd = (line_(m).mState) ...
                        && (line_(m).count < obj.totalNum) ...
                            && ~(line_(m).isBlock) && ~(line_(m).isStarve);

                        if isProd
                            line_(m).count = line_(m).count + 1;

                            if m == 1
                                line_(m).buffer_occ = line_(m).buffer_occ + 1;
                                CR(t) = CR(t) + 1;
                            elseif m == obj.M
                                line_(m - 1).buffer_occ = line_(m - 1).buffer_occ - 1;
                                PR(t) = PR(t) + 1;
                            else
                                line_(m).buffer_occ = line_(m).buffer_occ + 1;
                                line_(m - 1).buffer_occ = line_(m - 1).buffer_occ - 1;
                            end

                        end

                    end

                    % record WIP
                    for m = 1:obj.M - 1
                        WIP(m, t) = WIP(m, t) + line_(m).buffer_occ;
                    end

                    % record CT
                    c = cumsum(obj.B);

                    for k = 1:obj.batchNum

                        if line_(obj.M).count == c(k) && complete_flag(k) == 0
                            CT(k) = CT(k) + t;
                            complete_flag(k) = 1;
                        end

                    end

                end

            end

            % means
            PR = PR ./ ITER;
            CR = CR ./ ITER;
            WIP = WIP ./ ITER;
            ST = ST ./ ITER;
            BL = BL ./ ITER;
            CT = CT ./ ITER;
        end

        function [PR, CR, WIP, ST, BL, CT, T_RUN] = approximateAnalysis(obj, T, h_0, s_0, SingleMachineType, wip_cal_method)
            % aggregation based analysis method
            % T：evaluation step, if T==0, then automatically determine T and return it as T_RUN
            % h_0: initial buffer occupancy
            % s_0: initial machine status
            % SingleMachineType: The type of virtual one-machine line used for aggregation (Geo,Ber,Hyb)
            % wip_cal_method: WIP calculation method, 1=estimated from Bernoulli one-machine line 2=estimated from geometric one-machine line 3=estimated from geometric two-machine line
            if nargin < 5
                SingleMachineType = 'Geo';
            end

            if nargin < 6
                wip_cal_method = 1;
            end

            function [A, obj] = calStructTransMatrix(n)
                % claculate the template of A2
                doubleLine = GeometricDoubleMachineLine([5 7], [11 13], n);
                A = doubleLine.A;
                obj = doubleLine;
            end

            function A = calTransMatrixbyReplacement(As, p1, p2, r1, r2)
                % calculate A2 form template, i.e., through value replacement
                % AvalVirtual->Avalnew
                AvalVirtual = [24 -28 -30 35 -52 48 65 -60 -66 77 60 -70 143 -132 -130 120];
                Avalnew = [(p1 - 1) * (p2 - 1), p2 * (1 - p1), p1 * (1 - p2), p1 * p2, r2 * (1 - p1), (1 - p1) * (1 - r2), p1 * r2, p1 * (1 - r2), r1 * (1 - p2), p2 * r1, (1 - r1) * (1 - p2), p2 * (1 - r1), r1 * r2, r1 * (1 - r2), r2 * (1 - r1), (1 - r1) * (1 - r2)];
                A = zeros(size(As));

                for iAval = 1:16
                    A(As == AvalVirtual(iAval)) = Avalnew(iAval);
                end

            end

            function A = calGeoTransMatrix(p, r, B)
                % calculate the transition matrix of a geometric one-machine line
                S = 2 * (B + 1);
                A = zeros(S, S);

                for f = 0:B - 1
                    j = 2 * f + 1; %s=0
                    % repair
                    A(j + 1, j) = r;
                    % no repair
                    A(j, j) = 1 - r;

                    j = 2 * f + 2; %s=1
                    % breakdown
                    A(j + 1, j) = p;
                    % no breakdown
                    A(j + 2, j) = 1 - p;

                end

                A(2 * B + 1, 2 * B + 1) = 1;
                A(2 * B + 2, 2 * B + 2) = 1;
                A = sparse(A);

            end

            function A = calGeoTransMatrixbyReplacement(As, p, r)
                % A = zeros(size(As));
                A = As;
                A(As == 5) = r;
                A(As == -4) = 1 - r;
                A(As == 3) = p;
                A(As == -2) = 1 - p;
                % A(As == 1) = 1;
            end

            function A = calBerTransMatrix(p, B)
                % calculate the transition matrix of Bernoulli one-machine line
                A = eye(B + 1) * (1 - p);
                A = A + diag(ones(1, B), -1) * p;
                A(end, end) = 1;
                A = sparse(A);
            end

            function A = calBerTransMatrixbyReplacement(As, p)
                % A = zeros(size(As));
                A = As;
                A(As == 3) = p;
                A(As == -2) = 1 - p;
                % A(As == 1) = 1;
            end

            % pre-allocation of space
            T_PRE = 1000;
            PR = zeros(1, T_PRE);
            CR = zeros(1, T_PRE);
            WIP = zeros(obj.M - 1, T_PRE);
            ST = zeros(obj.M, T_PRE);
            BL = zeros(obj.M, T_PRE);
            CT = zeros(obj.batchNum, 1);

            % S.0 init
            % initialization of L2
            auxs = repmat(struct('A', [], 'X', []), obj.M - 1, 1);

            for m = 1:obj.M - 1
                [As, line] = calStructTransMatrix(obj.n(m));
                auxs(m).A = As;

                if nargin < 3
                    s1_init = 0;
                    s2_init = 0;
                    h_init = 0;
                else
                    isStarve = m ~= 1 && h_0(m - 1) == 0;
                    % recursive blocking is ignored, should not be a problem
                    isBlock = (m + 1) ~= obj.M && h_0(m) == obj.n(m) && s_0(m + 1) == 0;
                    s1_init = (s_0(m) && ~isStarve);
                    s2_init = (s_0(m + 1) && ~isBlock);
                    h_init = h_0(m);
                end

                idx_init = line.state2idx([h_init, s1_init, s2_init]);
                % fill in the initial state vector
                auxs(m).X = zeros(4 * (obj.n(m) + 1), 1);
                auxs(m).X(idx_init) = 1;
            end

            % initialization of L1
            % batch_auxs_geo = repmat(struct('A', [], 'C1', [], 'x', []), obj.M, 1);
            % batch_auxs_ber = repmat(struct('A', [], 'x', []), obj.M, 1);
            batch_auxs_geo = repmat(struct('A', calGeoTransMatrix(3, 5, obj.totalNum), ...
            'C1', [repmat([0 1], 1, obj.totalNum), 0, 0], ...
                'x', zeros(2 * (obj.totalNum + 1), 1)), obj.M, 1);
            batch_auxs_ber = repmat(struct('A', calBerTransMatrix(3, obj.totalNum), ...
                'x', zeros(obj.totalNum + 1, 1)), obj.M, 1);

            % Geo

            for i = 1:obj.M
                % batch_auxs_geo(i).x = zeros(2 * (obj.totalNum + 1), 1);
                % batch_auxs_geo(i).C1 = repmat([0 1], 1, obj.totalNum + 1);
                % batch_auxs_geo(i).C1(end) = 0;
                % batch_auxs_geo(i).A = calGeoTransMatrix(3, 5, obj.totalNum);

                if nargin > 2
                    % initialized by s_0 when initial state is provided
                    idx = 1 + 1 * (s_0(i) == 1);
                    batch_auxs_geo(i).x(idx) = 1;
                else
                    % set default machine state to 'off'
                    batch_auxs_geo(i).x(1) = 1;
                end

            end

            % Ber

            for i = 1:obj.M
                % batch_auxs_ber(i).x = zeros(obj.totalNum + 1, 1);
                % batch_auxs_ber(i).A = calBerTransMatrix(3, obj.totalNum);
                batch_auxs_ber(i).x(1) = 1;
            end

            % pre-allocation
            prs = zeros(1, obj.M - 1);
            crs = zeros(1, obj.M - 1);
            p_exp = zeros(1, obj.M);
            r_exp = zeros(1, obj.M);

            % pre-allocation
            ct_prob_sum = 0;
            t = 1;
            last_p_ct = 1;
            p_ct = 0;

            if T == 0
                % maximum allowed T
                MAX_T = 3000;
            else
                MAX_T = T;
            end

            while (abs(last_p_ct - p_ct) > 1e-6 || ct_prob_sum < 0.95) && t <= MAX_T
                % S.1 update p_exp, r_exp

                if strcmp(SingleMachineType, 'Geo')
                    % Geo mode
                    % update p_exp, r_exp from L1

                    for m = 1:obj.M
                        x_f = batch_auxs_geo(m).x(1:2:end) + batch_auxs_geo(m).x(2:2:end);
                        p_exp(m) = obj.efficienceMappingTable1(m, :) * x_f;
                        r_exp(m) = obj.efficienceMappingTable2(m, :) * x_f;
                    end

                else
                    % Hyb or Ber mode
                    % update p_exp, r_exp from L1

                    for m = 1:obj.M
                        p_exp(m) = obj.efficienceMappingTable1(m, :) * batch_auxs_ber(m).x;
                        r_exp(m) = obj.efficienceMappingTable2(m, :) * batch_auxs_ber(m).x;
                    end

                end

                pfs = p_exp;
                rfs = r_exp;
                pbs = p_exp;
                rbs = r_exp;

                % S.2 calculate pr, cr of L2
                for m = 1:obj.M - 1
                    C = repmat([0 1 0 1], 1, obj.n(m)); %Ip
                    C1 = [0 0 0 0 C];
                    C = repmat([0 0 1 1], 1, obj.n(m)); %Ic
                    C2 = [C 0 0 0 1];
                    prs(m) = C1 * auxs(m).X;
                    crs(m) = C2 * auxs(m).X;
                end

                % S.3 forward iteration
                for m = 2:obj.M
                    AA = calTransMatrixbyReplacement(auxs(m - 1).A, pfs(m - 1), p_exp(m), rfs(m - 1), r_exp(m));
                    C = repmat([0 1 0 1], 1, obj.n(m - 1)); %Ip
                    C1 = [0 0 0 0 C];
                    C2 = 1 - C1; %Inp
                    %pr=C1*auxs(m-1).X;
                    pr = prs(m - 1);

                    if pr == 0
                        pfs(m) = 1;
                    else
                        pfs(m) = (C2 * AA * (auxs(m - 1).X .* C1')) / pr;
                    end

                    if (1 - pr) == 0
                        rfs(m) = 1;
                    else
                        rfs(m) = (C1 * AA * (auxs(m - 1).X .* C2')) / (1 - pr);
                    end

                end

                %S.4 backward iteration
                for m = obj.M - 1:-1:1
                    AA = calTransMatrixbyReplacement(auxs(m).A, p_exp(m), pbs(m + 1), r_exp(m), rbs(m + 1));
                    C = repmat([0 0 1 1], 1, obj.n(m)); %Ic
                    C2 = [C 0 0 0 1];
                    C1 = 1 - C2; %Inc
                    %cr=C2*auxs(m).X;
                    cr = crs(m);

                    if cr == 0
                        pbs(m) = 0;
                    else
                        pbs(m) = (C1 * AA * (auxs(m).X .* C2')) / cr;
                    end

                    if (1 - cr) == 0
                        rbs(m) = 0;
                    else
                        rbs(m) = (C2 * AA * (auxs(m).X .* C1')) / (1 - cr);
                    end

                end

                %S.5 update p_hat, r_hat
                % use crs for the first machine and prs for the rest
                p_hat = zeros(1, obj.M);
                p_hat(2:end) = prs;
                p_hat(1) = crs(1);

                %S.6 update status vectors
                % for L1
                % Geo

                for m = 1:obj.M

                    if m == 1
                        A_Geo = calGeoTransMatrixbyReplacement(batch_auxs_geo(m).A, pbs(m), rbs(m));
                    else
                        A_Geo = calGeoTransMatrixbyReplacement(batch_auxs_geo(m).A, pfs(m), rfs(m));
                    end

                    batch_auxs_geo(m).x = A_Geo * batch_auxs_geo(m).x;

                end

                % Ber

                for m = 1:obj.M

                    A_Ber = calBerTransMatrixbyReplacement(batch_auxs_ber(m).A, p_hat(m));
                    batch_auxs_ber(m).x = A_Ber * batch_auxs_ber(m).x;

                end

                if strcmp(SingleMachineType, 'Geo') || strcmp(SingleMachineType, 'Hyb')
                    CR(t) = batch_auxs_geo(1).C1 * batch_auxs_geo(1).x;
                    PR(t) = batch_auxs_geo(obj.M).C1 * batch_auxs_geo(obj.M).x;
                else
                    CR(t) = [ones(1, obj.totalNum) * p_hat(1), 0] * batch_auxs_ber(1).x;
                    PR(t) = [ones(1, obj.totalNum) * p_hat(obj.M), 0] * batch_auxs_ber(obj.M).x;
                end

                % For L2
                for m = 1:obj.M - 1

                    AA = calTransMatrixbyReplacement(auxs(m).A, pfs(m), pbs(m + 1), rfs(m), rbs(m + 1));

                    nState = 4 * (obj.n(m) + 1);
                    C4 = zeros(1, nState); %st
                    C4(2) = 1; C4(4) = 1;
                    C5 = zeros(1, nState); %bl
                    C5(end - 1) = 1;

                    if strcmp(SingleMachineType, 'Geo')
                        probNotComplete = 1 - (batch_auxs_geo(m).x(end - 1) ...
                            + batch_auxs_geo(m).x(end)); % probability that a batch has not been compleated
                    else
                        probNotComplete = 1 - batch_auxs_ber(m).x(end);
                    end

                    switch wip_cal_method
                        case 1
                            % use Ber L1
                            % f(m)-f(m+1)
                            f1 = batch_auxs_ber(m).x;
                            f2 = batch_auxs_ber(m + 1).x;
                            WIP(m, t) = (0:length(f1) - 1) * (f1 - f2);
                            % WIP(m, t) = (1:length(batch_auxs_ber(m).x)) * batch_auxs_ber(m).x ...
                            %     - (1:length(batch_auxs_ber(m + 1).x)) * batch_auxs_ber(m + 1).x;
                        case 2
                            % use Geo L1
                            % f(m)-f(m+1)
                            f1 = batch_auxs_geo(m).x(1:2:end) + batch_auxs_geo(m).x(2:2:end);
                            f2 = batch_auxs_geo(m + 1).x(1:2:end) + batch_auxs_geo(m + 1).x(2:2:end);
                            WIP(m, t) = (0:length(f1) - 1) * (f1 - f2);
                        case 3
                            % use Geo L2
                            C3 = repmat((0:obj.n(m)), 4, 1); %wip
                            C3 = reshape(C3, [1, nState]);
                            WIP(m, t) = C3 * auxs(m).X * probNotComplete;
                        otherwise
                            error('Invalid wip_cal_method: %d', wip_cal_method);
                    end

                    ST(m + 1, t) = C4 * auxs(m).X * probNotComplete;
                    BL(m, t) = C5 * auxs(m).X * probNotComplete;
                    auxs(m).X = AA * auxs(m).X;
                end

                % claculate CT

                if T == 0
                    last_p_ct = p_ct;
                else
                    last_p_ct = -1; % when T is specified, the error condition is made never to be satisfied to avoid early exit
                end

                % completion probability at t
                if strcmp(SingleMachineType, 'Geo')
                    % Geo
                    % CT_i
                    batchNumCumSum = cumsum(obj.B);

                    for k = 1:obj.batchNum
                        C_ct_k = zeros(1, 2 * (obj.totalNum + 1));
                        C_ct_k(2 * batchNumCumSum(k)) = 1;
                        p_ct = C_ct_k * batch_auxs_geo(obj.M).x;

                        if k == obj.batchNum
                            ct_prob_sum = ct_prob_sum + p_ct;
                        end

                        CT(k) = CT(k) + t * p_ct;
                    end

                    t = t + 1;

                else
                    % SingleMachineType = Ber or Hyb
                    % CT_i
                    batchNumCumSum = cumsum(obj.B);

                    for k = 1:obj.batchNum
                        C_ct_k = zeros(1, obj.totalNum + 1);
                        C_ct_k(batchNumCumSum(k)) = p_hat(obj.M);
                        p_ct = C_ct_k * batch_auxs_ber(obj.M).x;

                        if k == obj.batchNum
                            ct_prob_sum = ct_prob_sum + p_ct;
                        end

                        CT(k) = CT(k) + t * p_ct;
                    end

                    t = t + 1;

                end

            end

            T_RUN = t - 1;
            PR = PR(1:T_RUN);
            CR = CR(1:T_RUN);
            WIP = WIP(:, 1:T_RUN);
            BL = BL(:, 1:T_RUN);
            ST = ST(:, 1:T_RUN);

        end

        function [CT, PRss] = calApprCTUsingSteadyPR(obj)
            % estimate CT using PRss
            PRss = zeros(1, obj.batchNum);

            function Q = calQ(P1, R1, P2, R2, N)
                % the Q function of geometric two-machine line
                % see paper "Identifying Bottlenecks in Serial Production Lines with Geometric Machines: Indicators and Rules"
                if N == 1
                    beta_2 = R1 + R2 - R1 * R2 - P2 * R1;
                    Q = P1 * beta_2 / ((R1 + R2 - R1 * R2) * (R1 + P1));
                else
                    alpha_1 = P1 + P2 - P1 * P2 - R1 * P2;
                    alpha_2 = P1 + P2 - P1 * P2 - R2 * P1;
                    beta_1 = R1 + R2 - R1 * R2 - P1 * R2;
                    beta_2 = R1 + R2 - R1 * R2 - P2 * R1;
                    sig = (alpha_2 * beta_1) / (alpha_1 * beta_2);
                    tempA = P1 * R2 * alpha_1 * alpha_2 * beta_2 * (P2 + beta_2);
                    tempB = P1 * R1 * R2 * alpha_2 * (beta_2^2 + P2 * (alpha_1 + beta_1) * (alpha_2 + 2 * beta_2));
                    tempC = 0;

                    for k = 2:N - 1
                        tempC = tempC + P1 * P2 * R1 * R2 * (alpha_2 + beta_2)^3 * sig^(k - 1);
                    end

                    tempD = P2 * R1 * alpha_1 * beta_2 * (R2 * (alpha_1 + beta_1) + alpha_2 * (P1 + R1)) * sig^(N - 1);
                    Q = (P1 * alpha_1 * alpha_2 * beta_2^2 * (R2 + P2)) / (tempA + tempB + tempC + tempD);
                end

            end

            function PR = calApprGeoSteadyPR(p, r, n)
                % calculate steady-state PR for geometric serial lines by steady-state aggregation
                % see paper "Identifying Bottlenecks in Serial Production Lines with Geometric Machines: Indicators and Rules"
                pfs = p;
                pbs = p;
                rfs = r;
                rbs = r;
                Mmachine = length(p);

                for k = 1:20 % 20 iterations are enough for convergence
                    % backward
                    for mm = Mmachine - 1:-1:1
                        pbs(mm) = p(mm) + r(mm) * calQ(pbs(mm + 1), rbs(mm + 1), pfs(mm), rfs(mm), n(mm));
                        rbs(mm) = r(mm) - r(mm) * calQ(pbs(mm + 1), rbs(mm + 1), pfs(mm), rfs(mm), n(mm));
                    end

                    % forward
                    for mm = 2:Mmachine
                        rfs(mm) = r(mm) - r(mm) * calQ(pfs(mm - 1), rfs(mm - 1), pbs(mm), rbs(mm), n(mm - 1));
                        pfs(mm) = p(mm) + r(mm) * calQ(pfs(mm - 1), rfs(mm - 1), pbs(mm), rbs(mm), n(mm - 1));
                    end

                end

                PR = rfs(end) / (pfs(end) + rfs(end));

            end

            CT = 0;

            for iBatch = 1:obj.batchNum

                PR_batch = calApprGeoSteadyPR(obj.p(iBatch, :), obj.r(iBatch, :), obj.n);

                PRss(iBatch) = PR_batch;
                CT = CT + obj.B(iBatch) / PR_batch;
            end

        end

    end

end
