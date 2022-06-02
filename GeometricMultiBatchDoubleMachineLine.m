classdef GeometricMultiBatchDoubleMachineLine
    % GeometricMultiBatchDoubleMachineLine

    properties
        p % machine repair rate (matrix of size K*2)
        r % machine failure rate (matrix of size K*2)
        n % buffer size
        batchNum % number of batches
        totalNum % number of total parts
        B % batch size (vector of size 1*K)
        S % number of syatem states
        efficienceMappingTable1 % lookup table for p
        efficienceMappingTable2 % lookup table for r
        A % transition matrix
        C % observation matrix
        C_type % flag for part type
        Vh % the number of feasible h for each f
        Vh_cumsum % cumsum of Vh
    end

    methods

        function obj = GeometricMultiBatchDoubleMachineLine(p, r, n, B)
            % constructor
            %   p: machine repair rate (matrix of size K*2)
            %   r: machine failure rate (matrix of size K*2)
            %   n: buffer size
            %   B: batch size (vector of size 1*K)
            obj.p = p;
            obj.r = r;
            obj.n = n;
            obj.B = B;
            assert(size(p, 1) == size(r, 1));
            assert(size(p, 1) == length(B));
            assert(length(n) == 1);
            obj.batchNum = length(B);
            obj.totalNum = sum(B);
            % prepare lookup table for p, r
            obj.efficienceMappingTable1 = repelem(obj.p, obj.B, 1);
            obj.efficienceMappingTable1(end + 1, :) = obj.efficienceMappingTable1(end, :);
            % obj.efficienceMappingTable1(end + 1, :) = 1;
            obj.efficienceMappingTable2 = repelem(obj.r, obj.B, 1);
            obj.efficienceMappingTable2(end + 1, :) = obj.efficienceMappingTable2(end, :);

            if obj.totalNum > n
                obj.S = 4 * ((n + 1) * (obj.totalNum + 1 - n) + (1 + n) * n / 2);
                obj.Vh = [ones(1, obj.totalNum + 1 - n) * (n + 1), n:-1:1];
            else
                obj.S = 4 * (1 + obj.totalNum + 1) * (obj.totalNum + 1) / 2;
                obj.Vh = (obj.totalNum + 1):-1:1;
            end

            obj.Vh_cumsum = [0, cumsum(obj.Vh)];

            % obj.efficienceMappingTable2(end + 1, :) = 0;

            % for s = 1:obj.S
            %     state = obj.idx2state(s);
            %     idx = obj.state2idx(state);
            %     assert(s == idx);
            %     disp([s, state, idx]);
            % end

            obj.A = obj.calTransMatrix();
            obj.C = zeros(6, obj.S); %PR CR Pct ST BL WIP

            for s = 1:obj.S
                state = obj.idx2state(s); %(s1,s2,h,f)
                s1 = state(1);
                s2 = state(2);
                h = state(3);
                f = state(4);
                obj.C(1, s) = s2 && h ~= 0 && f < obj.totalNum; % PR
                obj.C(2, s) = s1 && ~(h == obj.n && s2 == 0) && ((f + h) < obj.totalNum); %CR
                obj.C(3, s) = (f == (obj.totalNum - 1)) && s2 && h ~= 0; %Pct
                obj.C(4, s) = (f < obj.totalNum) && s2 && h == 0; %ST
                obj.C(5, s) = (f < obj.totalNum) && s1 && (h == obj.n && s2 == 0); %BL
                obj.C(6, s) = (f < obj.totalNum) * h; % WIP
            end

            obj.C_type = zeros(obj.batchNum, obj.S); % type k
            % segment corresponding to each batch, e.g. batchSeq(k)~batchSeq(k+1) for type k
            batchSeq = [0, cumsum(obj.B) - 1];

            for k = 1:obj.batchNum
                C_k = zeros(1, obj.S);

                for f_k = batchSeq(k):batchSeq(k + 1)
                    C_k((4 * f_k + 1):(4 * (obj.totalNum + 1)):end) = 1;
                    C_k((4 * f_k + 2):(4 * (obj.totalNum + 1)):end) = 1;
                    C_k((4 * f_k + 3):(4 * (obj.totalNum + 1)):end) = 1;
                    C_k((4 * f_k + 4):(4 * (obj.totalNum + 1)):end) = 1;
                end

                obj.C_type(k, :) = C_k;

            end

        end

        function A = calTransMatrix(obj)
            % Calculate the one-step transition matrix A
            A = sparse(obj.S, obj.S);

            if isa(obj.p(1), 'sym')
                A = sym(A); % when the input is a sym variable, the A matrix should be transformed into a sym matrix
            end

            for j = 1:obj.S
                state = obj.idx2state(j);
                s1 = state(1);
                s2 = state(2);
                h = state(3);
                f = state(4);

                if f ~= obj.totalNum % non-absorbing state
                    % mStates = dec2bin(2^2 - 1:-1:0) - '0';
                    mStates = [1 1; 1 0; 0 1; 0 0];

                    for imStates = 1:4
                        s1_new = mStates(imStates, 1);
                        s2_new = mStates(imStates, 2);
                        p_1 = obj.efficienceMappingTable1(f + h + 1, 1);
                        r_1 = obj.efficienceMappingTable2(f + h + 1, 1);
                        p_2 = obj.efficienceMappingTable1(f + 1, 2);
                        r_2 = obj.efficienceMappingTable2(f + 1, 2);
                        prob1 = (s1 == 1 && s1_new == 0) * p_1 + ...
                            (s1 == 1 && s1_new == 1) * (1 - p_1) + ...
                            (s1 == 0 && s1_new == 1) * r_1 + ...
                            (s1 == 0 && s1_new == 0) * (1 - r_1);
                        prob2 = (s2 == 1 && s2_new == 0) * p_2 + ...
                            (s2 == 1 && s2_new == 1) * (1 - p_2) + ...
                            (s2 == 0 && s2_new == 1) * r_2 + ...
                            (s2 == 0 && s2_new == 0) * (1 - r_2);
                        m2_prod = s2 && (f < obj.totalNum);
                        m1_prod = s1 && ((f + h) < obj.totalNum);
                        h_new = h - m2_prod * min(h, 1);
                        h_new = h_new + m1_prod * min(obj.n - h_new, 1);
                        f_new = f + ((h ~= 0) && m2_prod);

                        i = obj.state2idx([s1_new, s2_new, h_new, f_new]);
                        A(i, j) = prob1 * prob2;
                    end

                else
                    A(j, j) = 1;
                end

            end

        end

        function idx = state2idx(obj, state)
            % state=(s1,s2,h,f) -> index
            idx = 4 * obj.Vh_cumsum(state(4) + 1) + 4 * state(3) + 2 * state(1) + state(2) + 1;
        end

        function state = idx2state(obj, idx)
            % index -> state=(s1,s2,h,f)
            nh = floor((idx - 1) / 4) + 1;
            f = find(obj.Vh_cumsum < nh, 1, 'last') - 1;
            h = nh - obj.Vh_cumsum(f + 1) - 1;
            s1 = floor(mod((idx - 1), 4) / 2);
            s2 = mod(idx - 1, 2);
            state = [s1, s2, h, f];
        end

        function [PR, CR, ST, BL, CT, WIP] = markovAnalysis(obj, SIM_T, x_0)
            % Exact analysis with Markov analysis
            % SIM_T: number of evaluation steps
            % x_0 initial state (s,f)

            % set initial system state
            x = zeros(obj.S, 1);

            if nargin < 3
                % default: s1=0 s2=0 h=0 f=0
                init_idx = obj.state2idx([0, 0, 0, 0]);
            else
                init_idx = obj.state2idx(x_0);
            end

            x(init_idx) = 1;
            % initialization
            T_PRE = 3000; % pre-allocation size

            PR = zeros(1, T_PRE);
            CR = zeros(1, T_PRE);
            ST = zeros(1, T_PRE);
            BL = zeros(1, T_PRE);
            WIP = zeros(1, T_PRE);
            CT = 0;

            if SIM_T == 0
                MAX_T = T_PRE; % maximum allowed T
            else
                MAX_T = SIM_T;
            end

            t = 1;
            ct_prob_sum = 0;

            while (ct_prob_sum < 0.99999 || SIM_T ~= 0) && t <= MAX_T
                PR(t) = obj.C(1, :) * x;
                CR(t) = obj.C(2, :) * x;
                Pct = obj.C(3, :) * x;
                ct_prob_sum = ct_prob_sum + Pct;
                CT = CT + t * Pct;
                ST(t) = obj.C(4, :) * x;
                BL(t) = obj.C(5, :) * x;
                WIP(t) = obj.C(6, :) * x;
                x = obj.A * x;
                t = t + 1;
            end

            t = t - 1;
            PR = PR(1:t);
            CR = CR(1:t);
            ST = ST(1:t);
            BL = BL(1:t);
            WIP = WIP(1:t);
            % disp(t);

        end

        function [CT, Ect, Dct, Pabsorb] = calCT(obj, x_0)
            % ref: Introduction to Stochastic Processes with R, P119

            if nargin < 2
                x_0 = obj.state2idx([0, 0, 0, 0]);
            else
                x_0 = obj.state2idx(x_0);
            end

            AA = obj.A';
            Q = AA(1:end - 4, 1:end - 4); % 4 transient states
            R = AA(1:end - 4, end - 3:end);
            F = inv(eye(obj.S - 4) - Q); % fundamental matrix
            Pabsorb = F * R; % absorption probability
            Ect = F * ones(obj.S - 4, 1); %  mean absorption time of entering a certain absorption state
            Dct = (2 * F - eye(obj.S - 4)) * Ect - Ect.^2; % var
            CT = Ect(x_0);
        end

        function [CT, PRss] = calApprCTUsingSteadyPR(obj)
            % estimate CT using PRss
            % return CT and PRss of each batch (vector of size 1*K)

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

                % exact analysis
                PR = r(2) / (p(2) + r(2)) * (1 - calQ(p(1), r(1), p(2), r(2), n));
            end

            CT = 0;

            for iBatch = 1:obj.batchNum
                PR_batch = calApprGeoSteadyPR(obj.p(iBatch, :), obj.r(iBatch, :), obj.n);
                PRss(iBatch) = PR_batch;
                CT = CT + obj.B(iBatch) / PR_batch;
            end

        end

        function mc = showGraphplot(obj)
            % show state transfer diagrams
            % require Econometrics Toolbox
            mc = dtmc(obj.A');
            figure;
            graphplot(mc, 'ColorNodes', true);

        end

        function [PR, CR, CT] = approximateAnalysis(obj, SIM_T, x_0)
            % aggregation based analysis method
            % SIM_T: number of evaluation steps
            % x_0 initial state (s,f)

            % set initial system state
            x = zeros(obj.S, 1);

            if nargin < 3
                % default: s1=0 s2=0 h=0 f=0
                init_idx = obj.state2idx([0, 0, 0, 0]);
            else
                init_idx = obj.state2idx(x_0);
            end

            x(init_idx) = 1;
            % state vector of L1
            X = [1; 0];

            PR = zeros(1, SIM_T);
            CR = zeros(1, SIM_T);
            CT = 0;

            for t = 1:SIM_T
                pr = obj.C(1, :) * x;
                % cr = obj.C(2, :) * x;
                if pr == 0
                    p_f = 1;
                else
                    p_f = (1 - obj.C(1, :)) * obj.A * (x .* obj.C(1, :)') / pr;
                end

                if (1 - pr) == 0
                    r_f = 1;
                else
                    r_f = obj.C(1, :) * obj.A * (x .* (1 - obj.C(1, :))') / (1 - pr);
                end

                % calculate A_f and PR of L1
                A_f = [1 - r_f p_f; r_f 1 - p_f];
                PR(t) = X(2);
                X = A_f * X;

                % PR(t) = l.C1 * X;
                x = obj.A * x;
            end

        end

        function [PR, CR, CT] = MonteCarloAnalysis(obj, T, x_0, ITER)
            % Monte Carlo analysis of geometric lines
            % x_0 initial state (s,f)
            % ITER: number of repetitions

            if nargin < 4
                ITER = 1000;
            end

            PR = zeros(1, T);
            CR = zeros(1, T);
            CT = 0;

            for iITER = 1:ITER

                if nargin < 3
                    s = [0, 0];
                    h = 0;
                    f = 0;
                else
                    s = x_0(1:2);
                    h = x_0(3);
                    f = x_0(4);
                end

                complete_flag = 0;
                prod_flag = [0 0];

                for t = 1:T
                    p_rand = rand(1, 2);
                    r_rand = rand(1, 2);
                    p_1 = obj.efficienceMappingTable1(f + h + 1, 1);
                    r_1 = obj.efficienceMappingTable2(f + h + 1, 1);
                    p_2 = obj.efficienceMappingTable1(f + 1, 2);
                    r_2 = obj.efficienceMappingTable2(f + 1, 2);
                    p_ = [p_1, p_2];
                    r_ = [r_1, r_2];

                    prod_flag(2) = s(2) && f < obj.totalNum && h ~= 0;
                    prod_flag(1) = s(1) && (f + h) < obj.totalNum && ~(h == obj.n && ~s(2));

                    s = ((p_rand > p_) & (s == 1)) ...
                        | ((r_rand < r_) & (s == 0));

                    if prod_flag(2)
                        PR(t) = PR(t) + 1;
                        h = h - 1;
                        f = f + 1;
                    end

                    if prod_flag(1)
                        CR(t) = CR(t) + 1;
                        h = h + 1;
                    end

                    if f == obj.totalNum && complete_flag == 0
                        CT = CT + t;
                        complete_flag = 1;
                    end

                end

            end

            PR = PR / ITER;
            CR = CR / ITER;
            CT = CT / ITER;

        end

    end

end
