function max_yCI95 = plot_ci(y)
    % Plot the confidence interval of y, y is a k*N matrix, k is the number of experiments, N is the length of the data series
    % e.g.
    % x = 1:100;
    % y = sin(0.05*x)+0.5*cos(0.05*x).*randn(50,100);
    % plot_ci(y)

    N = size(y, 1);
    x = 1:size(y, 2);
    yMean = mean(y);
    ySEM = std(y) / sqrt(N);

    % since the overall variance was unknown, a t-test was used to calculate 95% CI intervals
    CI95 = tinv([0.025 0.975], N - 1);
    yCI95 = bsxfun(@times, ySEM, CI95(:));
    max_yCI95 = max(yCI95(2, :));

    % plot
    figure
    hold on
    plot(x, yMean, 'blue', 'LineWidth', 1.2)
    X_plot = [x, fliplr(x)];
    Y_plot = [yMean + yCI95(1, :), fliplr(yMean + yCI95(2, :))];
    fill(X_plot, Y_plot, 1, ...
        'facecolor', 'blue', ...
        'edgecolor', 'none', ...
        'facealpha', 0.3);
    grid
    % hold off
end
