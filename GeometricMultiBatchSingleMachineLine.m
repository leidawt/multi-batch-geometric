classdef GeometricMultiBatchSingleMachineLine
    % GeometricMultiBatchSingleMachineLine

    properties
        p % machine repair rate (matrix of size K*1)
        r % machine failure rate (matrix of size K*1)
        A % transition matrix
        batchNum % number of batches
        totalNum % number of total parts
        B % batch size
        S % number of syatem states
        efficienceMappingTable % lookup table
        C1 % observation matrix for PR
        C2 % observation matrix for Pct

    end

    methods

        function obj = GeometricMultiBatchSingleMachineLine(p, r, B)
            % constructor
            %   p: machine repair rate (matrix of size K*1)
            %   r: machine failure rate (matrix of size K*1)
            %   B: batch size
            obj.p = p;
            obj.r = r;
            obj.B = B;
            assert(size(p, 1) == size(r, 1));
            assert(size(p, 1) == length(B));
            obj.batchNum = length(B);
            obj.totalNum = sum(B);
            obj.S = 2 * (obj.totalNum + 1);
            % prepare lookup table
            obj.efficienceMappingTable = repelem([obj.p, obj.r], obj.B, 1);
            obj.efficienceMappingTable(end + 1, :) = obj.efficienceMappingTable(end, :);
            obj.efficienceMappingTable = obj.efficienceMappingTable';
            % prepare observation matrix for PR, Pct
            obj.C1 = repmat([0 1], 1, obj.totalNum + 1);
            obj.C1(end) = 0;
            obj.C2 = zeros(1, obj.S);
            obj.C2(end - 2) = 1;
            obj.A = obj.calTransMatrix();
        end

        function A = calTransMatrix(obj)
            % Calculate the one-step transition matrix A
            A = zeros(obj.S, obj.S);

            if isa(obj.p(1, 1), 'sym')
                A = sym(A); % when the input is a sym variable, the A matrix should be transformed into a sym matrix
            end

            for f = 0:obj.totalNum - 1
                p_ = obj.efficienceMappingTable(1, f + 1);
                r_ = obj.efficienceMappingTable(2, f + 1);
                j = 2 * f + 1; %s=0
                % repair
                A(j + 1, j) = r_;
                % non-repair
                A(j, j) = 1 - r_;

                j = 2 * f + 2; %s=1
                % breakdown
                A(j + 1, j) = p_;
                % non-breakdown
                A(j + 2, j) = 1 - p_;

            end

            A(2 * obj.totalNum + 1, 2 * obj.totalNum + 1) = 1;
            A(2 * obj.totalNum + 2, 2 * obj.totalNum + 2) = 1;

        end

        function [PR, CR, CT] = MonteCarloAnalysis(obj, T, x_0, ITER)
            % Monte Carlo analysis of geometric lines
            % x_0: initial state (s,f)
            % ITER: number of repetitions

            if nargin < 4
                ITER = 10000;
            end

            PR = zeros(1, T);
            CT = 0;

            for iITER = 1:ITER

                if nargin < 3
                    s = 0;
                    f = 0;
                else
                    s = x_0(1);
                    f = x_0(2);
                end

                complete_flag = 0;

                for t = 1:T
                    p_rand = rand(1);
                    r_rand = rand(1);
                    p_ = obj.efficienceMappingTable(1, f + 1);
                    r_ = obj.efficienceMappingTable(2, f + 1);

                    if s == 1 && f < obj.totalNum
                        PR(t) = PR(t) + 1;
                        f = f + 1;
                    end

                    s = (p_rand > p_ && s == 1) ...
                        || (r_rand < r_ && s == 0);

                    if f == obj.totalNum && complete_flag == 0
                        CT = CT + t;
                        complete_flag = 1;
                    end

                end

            end

            PR = PR / ITER;
            CR = PR;
            CT = CT / ITER;

        end

        function [CT, Ect, Dct, Pabsorb] = calCT(obj, x_0)
            % ref: Introduction to Stochastic Processes with R, P119
            % x_0: initial state (s,f)

            if nargin < 2
                x_0 = 1;
            else
                x_0 = x_0(1) + 2 * x_0(2) + 1;
            end

            AA = obj.A';
            Q = AA(1:end - 2, 1:end - 2); % 2 transient states
            R = AA(1:end - 2, end - 1:end);
            F = inv(eye(obj.S - 2) - Q); % fundamental matrix
            Pabsorb = F * R; % absorption probability
            Ect = F * ones(obj.S - 2, 1); % mean absorption time of entering a certain absorption state
            Dct = (2 * F - eye(obj.S - 2)) * Ect - Ect.^2; % var
            CT = Ect(x_0);
        end

        function CT = calCTFast(obj, x_0)
            % calculate CT by analytical formula
            % x_0: initial state (s,f)

            if nargin < 2
                x_0 = 1;
            end

            x_0 = x_0(1);
            CT = zeros(1, obj.batchNum);
            CT(1) = (obj.B(1) - 1) * obj.p(1) / obj.r(1) + obj.B(1) + (1 - x_0) / obj.r(1);

            for i = 2:obj.batchNum
                CT(i) = (obj.B(i) - 1) * obj.p(i) / obj.r(i) + obj.B(i) ...
                    + obj.p(i - 1) / obj.r(i);
            end

            CT = cumsum(CT);

        end

        function mc = showGraphplot(obj)
            % show state transfer diagrams
            % require Econometrics Toolbox
            mc = dtmc(obj.A');
            figure;
            graphplot(mc, 'ColorNodes', true);

        end

        function [PR, CR, CT] = markovAnalysis(obj, T, x_0)
            % Exact analysis with Markov analysis
            % T: number of evaluation steps
            % x_0: initial state (s,f)
            x = zeros(obj.S, 1);

            if nargin < 3
                x(1) = 1;
            else
                x(x_0(1) + 2 * x_0(2) + 1) = 1;
            end

            x_ = x;

            PR = zeros(1, T);

            for t = 1:T
                PR(t) = obj.C1 * x;
                x = obj.A * x;
            end

            CR = PR;

            probSum = 0;
            t = 1;
            CT = 1;

            while probSum < 0.9999999
                x_ = obj.A * x_;
                p_ct = obj.C2 * x_;
                CT = CT + t * p_ct;
                t = t + 1;
                probSum = probSum + p_ct(end);
            end

        end

    end

end
