classdef GeometricDoubleMachineLine
    % GeometricDoubleMachineLine

    properties
        p % machine repair rate (vector of size 1*2)
        r % machine failure rate (vector of size 1*2)
        n % buffer size
        S % number of syatem states
        A % transition matrix
        X % state vector
        hFirst
    end

    methods

        function obj = GeometricDoubleMachineLine(p, r, n, hFirst)
            % constructor
            %   p: machine repair rate (vector of size 1*2)
            %   r: machine failure rate (vector of size 1*2)
            %   n: buffer size
            assert(length(p) == 2);
            assert(length(r) == 2);
            assert(length(n) == 1);
            obj.p = p; obj.n = n; obj.r = r;
            obj.S = 4 * (n + 1);
            obj.X = zeros(obj.S, 1);

            if nargin < 4
                obj.hFirst = true; % first update h then s
            else
                obj.hFirst = hFirst;
            end

            obj.A = obj.calTransMatrix();

        end

        function A = calTransMatrix(obj)
            % Calculate the one-step transition matrix A
            A = sparse(obj.S, obj.S);

            if isa(obj.p(1), 'sym')
                A = sym(A); % when the input is a sym variable, the A matrix should be transformed into a sym matrix
            end

            for j = 1:obj.S
                state = obj.idx2state(j);
                % state = num2cell(state);
                % [h, s1, s2] = state{:};
                h = state(1); s1 = state(2); s2 = state(3);

                h_new = 0;

                if obj.hFirst
                    h_new = h - s2 * min(h, 1);
                    h_new = h_new + s1 * min(obj.n - h_new, 1);
                end

                % mStates = dec2bin(2^2 - 1:-1:0) - '0';
                mStates = [1 1; 1 0; 0 1; 0 0];

                for imStates = 1:length(mStates)
                    s1_new = mStates(imStates, 1);
                    s2_new = mStates(imStates, 2);
                    prob1 = (s1 == 1 && s1_new == 0) * obj.p(1) + (s1 == 1 && s1_new == 1) * (1 - obj.p(1)) + ...
                        (s1 == 0 && s1_new == 1) * obj.r(1) + (s1 == 0 && s1_new == 0) * (1 - obj.r(1));
                    prob2 = (s2 == 1 && s2_new == 0) * obj.p(2) + (s2 == 1 && s2_new == 1) * (1 - obj.p(2)) + ...
                        (s2 == 0 && s2_new == 1) * obj.r(2) + (s2 == 0 && s2_new == 0) * (1 - obj.r(2));

                    if ~obj.hFirst
                        h_new = h - s2_new * min(h, 1);
                        h_new = h_new + s1_new * min(obj.n - h_new, 1);
                    end

                    i = obj.state2idx([h_new, s1_new, s2_new]);
                    A(i, j) = prob1 * prob2;
                end

            end

            %obj.A=A;
        end

        function idx = state2idx(~, state)
            % state=(h, s1, s2) -> index
            idx = state(1) * 4 + state(2) * 2 + state(3) + 1;
        end

        function state = idx2state(~, idx)
            %index -> state=(h, s1, s2)
            h = floor((idx - 1) / 4);
            s1 = floor((idx - 1 - 4 * h) / 2);
            s2 = floor(idx - 1 - 4 * h - 2 * s1);
            state = [h, s1, s2];
        end

        function [pr, cr, wip, st, bl] = markovAnalysis(obj, SIM_T, x_0)
            % Exact analysis with Markov analysis
            % SIM_T: number of evaluation steps
            % x_0: initial state (s1, s2, h)
            C1 = repmat([0 1 0 1], 1, obj.n); %pr
            C1 = [0 0 0 0 C1];
            C2 = repmat([0 0 1 1], 1, obj.n); %cr
            C2 = [C2 0 0 0 1];
            C3 = repmat((0:obj.n), 4, 1); %wip
            C3 = reshape(C3, [1, obj.S]);
            C4 = zeros(1, obj.S); %st
            C4(2) = 1; C4(4) = 1;
            C5 = zeros(1, obj.S); %bl
            C5(end - 1) = 1;

            if nargin < 3
                x = zeros(obj.S, 1);
                init_idx = obj.state2idx([0, 0, 0]);
                x(init_idx) = 1; %default: s1=0 s2=0 h=0
            else
                x = x_0;
            end

            pr = zeros(SIM_T, 1);
            cr = zeros(SIM_T, 1);
            wip = zeros(SIM_T, 1);
            st = zeros(SIM_T, 1);
            bl = zeros(SIM_T, 1);

            for t = 1:SIM_T
                pr(t) = C1 * x;
                cr(t) = C2 * x;
                wip(t) = C3 * x;
                st(t) = C4 * x;
                bl(t) = C5 * x;
                x = obj.A * x;
            end

        end

    end

end
