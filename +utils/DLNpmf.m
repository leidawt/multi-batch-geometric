function [P, mean, std, cv] = DLNpmf(mu, sig)
    % calculate PMF for discrete log-normal distribution
    % https://doi.org/10.1007/s13198-021-01443-x
    % normcdf(x)
    x = 1:1000;
    P = normcdf((log(x) - mu) ./ sig) - normcdf((log(x - 1) - mu) ./ sig);
    mean = x * P';
    std = sqrt((x - mean).^2 * P');
    cv = std / mean;
    P = P(P > 1e-3);
    % disp([mean, std, cv]);
end
