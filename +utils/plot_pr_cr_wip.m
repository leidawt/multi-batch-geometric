function plot_pr_cr_wip(PR_appr, PR_sim, CR_appr, CR_sim, WIP_appr, WIP_sim, typeName)
    ax = nexttile;
    hold(ax, 'on')
    title(ax, sprintf('PR_{%s}', typeName));
    plot(ax, PR_appr, 'r');
    plot(ax, PR_sim, 'k');
    hold(ax, 'off');

    ax = nexttile;
    hold(ax, 'on')
    title(ax, sprintf('CR_{%s}', typeName));
    plot(ax, CR_appr, 'r');
    plot(ax, CR_sim, 'k');
    hold(ax, 'off');

    ax = nexttile;
    hold(ax, 'on')
    title(ax, sprintf('WIP_{%s}', typeName));
    plot(ax, WIP_appr(1, :), 'r');
    plot(ax, WIP_sim(1, :), 'k');
    hold(ax, 'off');
end
