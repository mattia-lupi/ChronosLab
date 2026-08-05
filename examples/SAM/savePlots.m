function savePlots(Name, sizeSeq, iter, t_setup, t_solve, params)
    % SAVEPLOTS Generates performance plot headlessly and exports image
    outputDir = fullfile('images', Name);
    if ~exist(outputDir, 'dir')
       mkdir(outputDir);
    end

    fig = figure('Name', "SAM Benchmark - " + Name, 'Color', 'w', ...
                 'Position', [100, 100, 1200, 800], 'Visible', 'off');

    labels = {'Recompute Precond', 'Reuse Precond', 'SAM Adaptive'};
    colors = [0.8500 0.3250 0.0980; 
              0.9290 0.6940 0.1250; 
              0.4660 0.6740 0.1880];
    steps = 1:sizeSeq;
    tlayout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tlayout, "Performance Analysis: Dataset [" + Name + "]", 'FontSize', 14, 'FontWeight', 'bold');

    % Panel 1: GMRES Iterations
    nexttile;
    hold on; 
    plot(steps, iter(:, 1), '-o', 'LineWidth', 1.8, 'Color', colors(1, :), 'MarkerSize', 3, 'MarkerFaceColor', colors(1, :));
    plot(steps, iter(:, 2), '-o', 'LineWidth', 1.8, 'Color', colors(2, :), 'MarkerSize', 3, 'MarkerFaceColor', colors(2, :));
    plot(steps, iter(:, 3), '-o', 'LineWidth', 1.8, 'Color', colors(3, :), 'MarkerSize', 3, 'MarkerFaceColor', colors(3, :));
    if isfield(params, 'q') && ~isempty(params.q)
       xline(params.q, ':', 'Color', colors(2, :), 'LineWidth', 1.5);
       labels{end+1} = 'Recomp. Precond (Reuse)';
    end
    if isfield(params, 'tP') && ~isempty(params.tP)
       xline(params.tP, ':', 'Color', colors(3, :), 'LineWidth', 1.5);
       labels{end+1} = 'Recomp. Precond (SAM)';
    end
    if isfield(params, 'tS') && ~isempty(params.tS)
       xline(params.tS, '--', 'Color', colors(3, :), 'LineWidth', 1.5);
       labels{end+1} = 'Compute SAM';
    end
    hold off; grid on;
    title('GMRES Iterations'); xlabel('Sequence Step'); ylabel('Iterations');
    legend(labels, 'Location', 'best');

    % Panel 2: Setup Time
    nexttile;
    b_setup = bar(steps, t_setup, 'grouped');
    for k = 1:3
       b_setup(k).FaceColor = colors(k, :);
    end
    grid on; title('Setup Time (Preconditioner + SAM)');
    xlabel('Sequence Step'); ylabel('Time (s)');
    legend(labels(1:3), 'Location', 'best');

    % Panel 3: Solve Time
    nexttile; hold on;
    for k = 1:3
       plot(steps, t_solve(:, k), '-s', 'LineWidth', 1.8, 'Color', colors(k, :), ...
            'MarkerSize', 3, 'MarkerFaceColor', colors(k, :));
    end
    hold off; grid on; % [cite: 52, 53]
    title('GMRES Solve Time'); xlabel('Sequence Step'); ylabel('Time (s)'); % [cite: 53]
    legend(labels, 'Location', 'best'); % [cite: 53]

    % Panel 4: Cumulative Total Time
    nexttile; start = 1; % [cite: 54]
    t_total = cumsum(t_setup + t_solve, 1); hold on; % [cite: 55]
    for k = start:3 % [cite: 55]
       plot(steps, t_total(:, k), '-^', 'LineWidth', 2.0, 'Color', colors(k, :), ...
            'MarkerSize', 3, 'MarkerFaceColor', colors(k, :)); % [cite: 55]
    end
    hold off; grid on; % [cite: 56]
    title('Cumulative Runtime (Setup + Solve)'); xlabel('Sequence Step'); % [cite: 56]
    ylabel('Total Accumulated Time (s)'); legend(labels(start:3), 'Location', 'northwest'); % [cite: 56]

    exportpath = fullfile(outputDir, "benchmark_results.png");
    saveas(fig, exportpath);
    close(fig);
    fprintf('Saved plot to %s\n', exportpath);
end