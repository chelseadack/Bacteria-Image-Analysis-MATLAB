function func_plot_smoothed_histogram_bacteria(MES, manualColors, scaleData, edges, smoothing, positive_threshold)
% FUNC_PLOT_SMOOTHED_HISTOGRAM_BACTERIA  Smoothed histogram of the fluorescence per cell, one line per condition, with the gate.
%
%   func_plot_smoothed_histogram_bacteria(MES, manualColors, scaleData, edges, smoothing, positive_threshold)
%
%   Plots into the current axes, so it can be used for a whole figure or for one tile of a tiledlayout.
%   The y axis is the % of cells in that condition in each bin.
%
%   Input parameters:
%       - MES = 1xC cell array, fluorescence of every cell in each condition.
%       - manualColors = 1xC cell array, one RGB colour per condition.
%       - scaleData = function handle used to scale the fluorescence values (eg. asinh).
%       - edges = histogram bin edges, on the scaled axis.
%       - smoothing = how much to smooth the lines, in bins (higher = smoother).
%       - positive_threshold = GFP positive gate, drawn as a dashed line.

    centres = edges(1:end-1) + diff(edges) / 2;

    hold on
    for j = 1:length(MES)
        counts = histcounts(scaleData(MES{j}), edges, 'Normalization', 'probability');
        plot(centres, 100 * smoothdata(counts, 'gaussian', smoothing), 'Color', manualColors{j}, 'LineWidth', 2.5)
    end
    xline(scaleData(positive_threshold), '--k', 'LineWidth', 1.5);

end
