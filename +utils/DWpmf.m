function [P, mean, std, cv] = DWpmf(p, beta)
    % calculate the PMF of the discrete Weibull distribution
    % http://www.math.wm.edu/~leemis/chart/UDR/PDFs/Discreteweibull.pdf
    % note: to match the shifted geometric distribution, a substitution, x=x-1, is made here
    x = 1:1000;
    P = (1 - p).^((x - 1).^beta) - (1 - p).^(x.^beta);
    mean = x * P';
    std = sqrt((x - mean).^2 * P');
    cv = std / mean;
    P = P(P > 1e-3);
    % disp([mean, std, cv]);
end
