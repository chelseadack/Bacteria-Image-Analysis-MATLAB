function func_plot_swarm_bacteria(MES, manualColors, scaleData, n, positive_threshold)
% FUNC_PLOT_SWARM_BACTERIA  Swarm plot of the fluorescence of every cell, one column per condition, with the gate.
%
%   func_plot_swarm_bacteria(MES, manualColors, scaleData, n, positive_threshold)
%
%   Plots into the current axes, so it can be used for a whole figure or for one tile of a tiledlayout.
%
%   Input parameters:
%       - MES = 1xC cell array, fluorescence of every cell in each condition.
%       - manualColors = 1xC cell array, one RGB colour per condition.
%       - scaleData = function handle used to scale the fluorescence values (eg. asinh).
%       - n = spacing between conditions on the x axis.
%       - positive_threshold = GFP positive gate, drawn as a dashed line.

    hold on
    for j = 1:length(MES)
        if ~isempty(MES{j}) % skip conditions where no cells were detected
            y = scaleData(MES{j}(:));

            % individual cells, spread horizontally
            swarmchart(j*n * ones(size(y)), y, 10, manualColors{j}, 'filled', 'MarkerFaceAlpha', 0.05, 'XJitter', 'density', 'XJitterWidth', 0.8)

            % black boxplot on top
            %boxchart(j*n * ones(size(y)), y, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0, 'WhiskerLineColor', 'k', ...
            %    'MarkerStyle', 'none', 'LineWidth', 1.5, 'BoxWidth', 0.6)
            %plot(j*n, mean(y), 'k*', 'MarkerSize', 10) % UNCOMMENT to mark the mean
        end
    end

    yline(scaleData(positive_threshold), '--k', 'LineWidth', 1.5);
    xlim([0 (length(MES)*n)+1])
    xticks(n * (1:length(MES)))

end
