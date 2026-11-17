%% CHELSEA'S BACTERIA ANALYSIS SCRIPT

% This script detects bacteria in brightfield microscope images and measures the average whole cell GFP fluorescence
% of every cell. GFP positive cells are found by gating against the negative control (Inactive Template). Outputs:
%	- per-cell fluorescence for each condition, for each repeat and for the whole population
%	  (swarm plots and histograms, with the GFP positive gate)
%	- % GFP positive cells and fold change vs the negative control, one data point per repeat
%	- switch response normalised to the negative (0%) and positive (100%) controls
% Various checkpoints are integrated which can be used to troubleshoot and check accuracy.

% Data must be in this format:
%	- One folder per condition (eg. '-AMF', '+AMF', 'Full Template', 'Inactive Template'),
%	  each containing all of the repeats for that condition.
%	- Colour or greyscale images must be saved as '.png'.
%	- Each unique image must have 2 images named 'repeat.image. condition':
%		Eg. 'B.3. -AMF_BF.png'  = Brightfield image used to detect the cells (repeat B, image 3).
%		    'B.3. -AMF_GFP.png' = Image of fluorescence to measure.
%	- The folder names will be used as graph labels as default, but these can be edited.

% To use this script:
%	1. Open in MATLAB online or in MATLAB on your PC (requires the Image Processing Toolbox).
%	2. Put the condition folders in the same folder as these scripts.
%	3. Run each section in order (Ctrl+Enter / Cmd+Enter), following the instructions in the code below.

%% ----------------- INPUT PARAMETERS -------------------------------------------------------------------------------------

% Close and erase current workspace, ensures only data it uses is from your current data.
close all
clc
clearvars

% define working directory (the folder containing these scripts and the condition folders) CHANGE THIS
working_folder = '/MATLAB Drive/MATLAB_CHELS_bacteria_quantification';
if ~isfolder(working_folder)
    error('working_folder not found: %s\nChange working_folder to the folder containing these scripts and your data.', working_folder);
end
addpath(working_folder)
addpath(fullfile(working_folder, 'Functions')) % folder containing the func_*.m files
cd(working_folder)

% define the condition folders to process, in the order they will be plotted
condition_folders = {'Inactive Template', 'Full Template', '-AMF', '+AMF'};

% Input the parameters used to detect the cells:
adapt_sensitivity = 5; % Noise factor: higher = fewer pixels called cell, lower = more pixels called cell
min_cell_area = 150; % Minimum cell area in pixels
max_cell_area = 2800; % Maximum cell area in pixels

% Input the parameters used to gate GFP positive cells:
negative_control = 'Inactive Template'; % condition used to set the gate
gate_percent = 5; % % of negative control cells allowed above the gate (5 or 1 are typical) CHANGE THIS

% Plot settings, shared by all of the figures below:
n = 1.5; % spacing between conditions on the x axis
jitter_width = 0.3; % horizontal spread of the individual data points CHANGE THIS

scale_factor = 1; % CHANGE THIS
scaleData = @(x) asinh(x / scale_factor); % asinh scale: linear near 0, log-like for large values

fluorescence_range = [0 255]; % fluorescence range to show on the swarm plots and histograms (8-bit images = 0 to 255) CHANGE THIS
tickValues = [0 1 2 3 5 10 20 50 100 200]; % fluorescence values to label on the asinh axis CHANGE THIS

% one colour per condition, in the same order as condition_folders
Grey = [124, 125, 107] / 255;
OriginBlue = [52, 100, 235] / 255;
OriginPurple = [163, 91, 212] / 255;
DarkGrey = [83, 82, 84] / 255;
manualColors = {Grey, OriginPurple, DarkGrey, OriginBlue};

% make axis labels the condition folder names
cleanNames = cellfun(@(x) strrep(x, '_', '\_'), condition_folders, 'UniformOutput', false); % escape underscores for plotting
%cleanNames = {'Inactive Template', 'Full Template', '- AMF', '+ AMF'}; % UNCOMMENT to manually define x axis labels

% NOW RUN THIS SECTION ^

%% ----------------- TEST PARAMETERS -------------------------------------------------------------------------------------
% Check the cell detection parameters are suitable for your samples and images.

% For a specified test sample, a figure will be show where it detected cells and where it will measure fluorescence.
% Specify the condition folder and the individual image you want to test (file name without '_BF.png'):
test_condition = '-AMF';
test_sample = 'B.1. -AMF';

% Create the figure
func_test_detection_bacteria(test_condition, test_sample, adapt_sensitivity, min_cell_area, max_cell_area)

% NOW RUN THIS SECTION ^
% Adjust parameters defined above accordingly.

%% ----------------- DETECT & MEASURE FLUORESCENCE FOR ALL SAMPLES -------------------------------------------------------------------------------------
% For all image files in every condition folder, detect cells and quantify fluorescence.

MES = cell(1, length(condition_folders)); % pooled values for each condition, all repeats
allMES = cell(1, length(condition_folders)); % values for each individual image in each condition
image_names = cell(1, length(condition_folders)); % image names in each condition
image_repeats = cell(1, length(condition_folders)); % repeat each image belongs to, in each condition

for f = 1:length(condition_folders)

    % Load in all image files for this condition
    [image_names{f}, image_repeats{f}, grey_DYE_images, grey_MES_images] = func_load_images_bacteria(condition_folders{f});

    % Detect cells and measure average fluorescence
    allMES{f} = func_detect_measure_bacteria(adapt_sensitivity, min_cell_area, max_cell_area, grey_DYE_images, grey_MES_images);

    % append data from every image in this condition together
    MES{f} = [allMES{f}{:}];

end

%clear grey_DYE_images grey_MES_images % free memory, the images are no longer needed

% pool the cells from every image in each repeat, for each condition
repeats = unique([image_repeats{:}]); % eg. {'A', 'B', 'C'}
MES_repeat = cell(1, numel(repeats)); % MES_repeat{r}{f} = values for repeat r, condition f
for r = 1:numel(repeats)
    MES_repeat{r} = cell(1, length(condition_folders));
    for f = 1:length(condition_folders)
        MES_repeat{r}{f} = [allMES{f}{strcmp(image_repeats{f}, repeats{r})}];
    end
end

% NOW RUN THIS SECTION ^
% The measured fluorescence values will be stored in 'MES' variable (one entry per condition, all repeats pooled)
% The values for each repeat will be stored in 'MES_repeat' variable (one entry per repeat, then one per condition)
% The values for each individual image will be stored in 'allMES' variable

%% ----------------- SET THE GFP POSITIVE GATE FROM THE NEGATIVE CONTROL -------------------------------------------------------------------------------------
% A cell is GFP positive if it is brighter than all but gate_percent % of the negative control cells.
% One gate is set from the negative control cells of all repeats pooled, and used for every repeat and every plot.

isNegative = strcmp(condition_folders, negative_control);
if ~any(isNegative)
    error('negative_control "%s" is not one of the condition_folders.', negative_control);
end

positive_threshold = 3;%func_gate_threshold(MES{isNegative}, gate_percent);
fprintf('GFP positive gate = %.2f A.U.\n', positive_threshold)

% % GFP positive cells in every image, and the average of the images in each repeat
[percentTable, repeatTable] = func_percent_positive_table(condition_folders, image_names, image_repeats, allMES, repeats, positive_threshold);

disp(percentTable)
disp(repeatTable)
% writetable(percentTable, 'percent_positive_per_image.csv'); % UNCOMMENT to save as a spreadsheet
% writetable(repeatTable, 'percent_positive_per_repeat.csv'); % UNCOMMENT to save as a spreadsheet

% NOW RUN THIS SECTION ^

%% ----------------- SWARM PLOT OF FLUORESCENCE PER CELL, WHOLE POPULATION -------------------------------------------------------------------------------------
% plot measured fluorescence values for every cell in each condition, all repeats pooled
%       - each dot is one cell, spread horizontally by local density
%       - dashed line = GFP positive gate
% edit y axis limits, font size, titles, colours, etc below as desired

figure2 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);
func_plot_swarm_bacteria(MES, manualColors, scaleData, n, positive_threshold)

ylabel({'Average Fluorescence Intensity', '(A.U., asinh scale)'}) % define y axis label
ylim(scaleData(fluorescence_range)) % COMMENT for automatic defined y axis limits
yticks(scaleData(tickValues)); yticklabels(string(tickValues));
xticklabels(cleanNames);
xtickangle(30)

ax1 = gca;
ax1.YAxis.FontSize = 18;
ax1.XAxis.FontSize = 18;
ax1.LineWidth=2;
grid off;
box off

%% ----------------- SWARM PLOT OF FLUORESCENCE PER CELL, EACH REPEAT -------------------------------------------------------------------------------------
% one panel per repeat, dashed line = GFP positive gate

figureSR = figure('Renderer', 'painters', 'Position', [10 10 400*numel(repeats) 700]);
t = tiledlayout(1, numel(repeats), 'TileSpacing', 'compact');
for r = 1:numel(repeats)
    nexttile
    func_plot_swarm_bacteria(MES_repeat{r}, manualColors, scaleData, n, positive_threshold)
    title(sprintf('Repeat %s', repeats{r}), 'FontWeight', 'normal')
    ylim(scaleData(fluorescence_range)) % COMMENT for automatic defined y axis limits
    yticks(scaleData(tickValues)); yticklabels(string(tickValues));
    if r > 1, yticklabels([]), end
    xticklabels(cleanNames); xtickangle(30)
    ax = gca; ax.FontSize = 16; ax.LineWidth = 2; box off
end
ylabel(t, {'Average Fluorescence Intensity', '(A.U., asinh scale)'}, 'FontSize', 18)

%% ----------------- HISTOGRAM OF FLUORESCENCE PER CELL, ONE PANEL PER CONDITION, WHOLE POPULATION --------DELETEEEEEEEEEE-------------------------------------------------------------------------------------

figureH = figure('Renderer', 'painters', 'Position', [10 10 500 800]);
t = tiledlayout(length(MES), 1, 'TileSpacing', 'compact');
for j = 1:length(MES)
    nexttile
    histogram(scaleData(MES{j}), 'BinWidth', 0.1, 'Normalization', 'probability', 'FaceColor', manualColors{j}, 'EdgeColor', 'none')
    ylim([0 0.03])
    xline(scaleData(positive_threshold), '--k', 'LineWidth', 1.5)
    text(0.98, 0.85, sprintf('%.0f%% positive', 100 * mean(MES{j} > positive_threshold)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 14)
    title(cleanNames{j}, 'FontWeight', 'normal')
    xticks(scaleData(tickValues))
    if j < length(MES), xticklabels([]), else, xticklabels(string(tickValues)), end
    ax = gca; ax.FontSize = 14; ax.LineWidth = 1.5; box off
end
linkaxes(findobj(t, 'Type', 'axes'), 'x') % same x axis for every panel
xlim(scaleData(fluorescence_range)) % fluorescence range to show
xlabel(t, {'Mean GFP fluorescence per cell', '(A.U., asinh scale)'}, 'FontSize', 16)
ylabel(t, 'Fraction of cells', 'FontSize', 16)

%% ----------------- SMOOTHED HISTOGRAMS OF FLUORESCENCE PER CELL, WHOLE POPULATION -------------------------------------------------------------------------------------
% all conditions on one graph, all repeats pooled, dashed line = GFP positive gate

binWidth = 0.05; % histogram bin width on the asinh scale
smoothing = 5; % how much to smooth the lines, in bins (higher = smoother) CHANGE THIS
edges = scaleData(fluorescence_range(1)) : binWidth : scaleData(fluorescence_range(2)) + binWidth;

figureS = figure('Renderer', 'painters', 'Position', [10 10 700 500]);
func_plot_smoothed_histogram_bacteria(MES, manualColors, scaleData, edges, smoothing, positive_threshold)

xticks(scaleData(tickValues)); xticklabels(string(tickValues));
xlim(scaleData(fluorescence_range)) % fluorescence range to show
ylim([0 1]) % y axis maximum CHANGE THIS (the negative peak will go off the top)
%set(gca, 'YScale', 'log') % UNCOMMENT to show the full negative peak and the tail together (then comment out ylim)
xlabel({'Mean GFP fluorescence per cell', '(A.U., asinh scale)'}); ylabel('Fraction of cells (%)');
legend(cleanNames, 'Location', 'northeast'); legend boxoff
ax = gca; ax.FontSize = 16; ax.LineWidth = 2; ax.Toolbar.Visible = 'off'; box off

%% ----------------- SMOOTHED HISTOGRAMS OF FLUORESCENCE PER CELL, EACH REPEAT -------------------------------------------------------------------------------------
% one panel per repeat, dashed line = GFP positive gate (uses binWidth, smoothing & edges from the section above)

figureSHR = figure('Renderer', 'painters', 'Position', [10 10 500*numel(repeats) 500]);
t = tiledlayout(1, numel(repeats), 'TileSpacing', 'compact');
for r = 1:numel(repeats)
    nexttile
    func_plot_smoothed_histogram_bacteria(MES_repeat{r}, manualColors, scaleData, edges, smoothing, positive_threshold)
    title(sprintf('Repeat %s', repeats{r}), 'FontWeight', 'normal')
    xticks(scaleData(tickValues)); xticklabels(string(tickValues));
    xlim(scaleData(fluorescence_range)) % fluorescence range to show
    ylim([0 1.5]) % y axis maximum CHANGE THIS (the negative peak will go off the top)
    if r > 1, yticklabels([]), end
    ax = gca; ax.FontSize = 16; ax.LineWidth = 2; ax.Toolbar.Visible = 'off'; box off
end
legend(cleanNames, 'Location', 'northeast'); legend boxoff % legend on the last panel only
xlabel(t, {'Mean GFP fluorescence per cell', '(A.U., asinh scale)'}, 'FontSize', 18)
ylabel(t, 'Fraction of cells (%)', 'FontSize', 18)

%% ----------------- SWITCH RESPONSE, NORMALISED TO THE CONTROLS -------------------------------------------------------------------------------------
% 0% = average of the negative control repeats (Inactive Template), 100% = average of the positive control repeats (Full Template)
% each dot = one repeat (the average of the images in that repeat), black line & error bar = mean ± SD of the repeats

show = condition_folders(ismember(condition_folders, {'-AMF', '+AMF'})); % conditions to plot, use show = condition_folders; to include the controls
%show = condition_folders;

pct = @(name) repeatTable.PercentPositive(strcmp(repeatTable.Condition, name));
negMean = mean(pct('Inactive Template'));
posMean = mean(pct('Full Template'));

figure6 = figure('Renderer', 'painters', 'Position', [10 10 400 600]);
hold on
for k = 1:numel(show)
    values = (pct(show{k}) - negMean) / (posMean - negMean) * 100;
    scatter(k + 0.2 * (rand(size(values)) - 0.5), values, 60, manualColors{strcmp(condition_folders, show{k})}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7)
    plot(k + [-0.25 0.25], mean(values) * [1 1], 'k', 'LineWidth', 2)
    errorbar(k, mean(values), std(values), 'k', 'LineWidth', 1.5, 'CapSize', 12)
end
%yline(0, ':k', 'Inactive Template', 'FontSize', 12); yline(100, ':k', 'Full Template', 'FontSize', 12);
yline(0, '--'); yline(100, '--');
xlim([0.5 numel(show) + 0.5]); xticks(1:numel(show)); xticklabels(strrep(show, '_', '\_'));
ylabel({'GFP positive cells', '(% of Full Template response)'});
ylim([-100 300])
ax = gca; ax.FontSize = 16; ax.LineWidth = 2; box off

%% ----------------- SWITCH RESPONSE, NORMALISED TO THE CONTROLS, INDIVIDUALLY -------------------------------------------------------------------------------------
% 0% = average of the negative control repeats (Inactive Template), 100% = average of the positive control repeats (Full Template)
% each dot = one repeat (the average of the images in that repeat), black line & error bar = mean ± SD of the repeats

repeatColors = [0.00 0.45 0.70;   % repeat A, blue
    0.84 0.37 0.00;   % repeat B, dark orange
    0.00 0.62 0.45];  % repeat C, green

show = condition_folders(ismember(condition_folders, {'-AMF', '+AMF'})); % conditions to plot, use show = condition_folders; to include the controls

% % GFP positive of one condition in one repeat (the average of that repeat's images)
pctRepeat = @(name, rep) repeatTable.PercentPositive(strcmp(repeatTable.Condition, name) & strcmp(repeatTable.Repeat, rep));

figure6 = figure('Renderer', 'painters', 'Position', [10 10 400 600]);
hold on
for k = 1:numel(show)
    % normalise each repeat to its own controls: 0% = that repeat's Inactive Template, 100% = that repeat's Full Template
    values = [];
    for r = 1:numel(repeats)
        neg = pctRepeat('Inactive Template', repeats{r});
        pos = pctRepeat('Full Template', repeats{r});
        x = pctRepeat(show{k}, repeats{r});
        if isempty(neg) || isempty(pos) || isempty(x), continue, end % skip repeats missing a control
        values(end+1) = (x - neg) / (pos - neg) * 100; %#ok<SAGROW>

        % one dot per repeat, coloured by repeat
        scatter(k + 0.2 * (rand - 0.5), values(end), 60, repeatColors(r,:), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('Repeat %s', repeats{r}))
    end
    plot(k + [-0.25 0.25], mean(values) * [1 1], 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
    errorbar(k, mean(values), std(values), 'k', 'LineWidth', 1.5, 'CapSize', 12, 'HandleVisibility', 'off')
end

%yline(0, ':k', 'Inactive Template', 'FontSize', 12); yline(100, ':k', 'Full Template', 'FontSize', 12);
yline(0, '--'); yline(100, '--');
[~, firstDots] = unique(get(findobj(gca, 'Type', 'scatter'), 'DisplayName'));
dots = findobj(gca, 'Type', 'scatter');
legend(dots(firstDots), 'Location', 'northwest'); legend boxoff
xlim([0.5 numel(show) + 0.5]); xticks(1:numel(show)); xticklabels(strrep(show, '_', '\_'));
ylabel({'GFP positive cells', '(% of Full Template response)'});
ylim([-50 250])
ax = gca; ax.FontSize = 16; ax.LineWidth = 2; box off
%% ----------------- OPTIONAL: BOX PLOT OF FLUORESCENCE PER CELL (func_al_goodplot, violins removed) -------------------------------------------------------------------------------------

figure7 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);

    for j = 1:length(MES)
        if ~isempty(MES{j}) % skip conditions where no cells were detected
            func_al_goodplot(scaleData(MES{j}),j*n, [], manualColors{j});
            %delete(findobj(gca, 'Type', 'scatter')) % hide individual data points
            set(findobj(gca, 'Type', 'scatter'), 'MarkerFaceAlpha', 0.05) % 0 = invisible, 1 = solid
        end
    end

    delete(findobj(gca, 'Type', 'patch', 'FaceAlpha', 0.15)) % remove the violins
    points = findobj(gca, 'Type', 'scatter');
    for k = 1:numel(points)
        points(k).XData = points(k).XData + jitter_width * (rand(size(points(k).XData)) - 0.5); % spread points horizontally
    end

    ylabel('Average Fluorescence Intensity (A.U.)') % define y axis label
    %ylim(scaleData(fluorescence_range)) % UNCOMMENT for manually defined y axis limits

    xlim([0 (length(MES)*n)+1])
    xticks(n * (1:length(MES)))
    xticklabels(cleanNames);
    xtickangle(30)

    ax1 = gca;
    ax1.YAxis.FontSize = 18;
    ax1.XAxis.FontSize = 18;
    ax1.YAxis.Exponent = 0;
    ax1.LineWidth=2;
    box off

% NOW RUN THIS SECTION ^
% Use the below section to check assigned values for individual samples & cells
% Adjust detection parameters, images, etc & repeat process, if desired
% Plots can be saved as pngs with the SAVE FIGURES section at the end of this script

%% ----------------- CHECK MEASURED FLUORESCENCE IS AS EXPECTED -------------------------------------------------------------------------------------

% This will overlay a coloured circle onto each cell in the specified image, indicating its measured fluorescence.
% Used to check the cells are being quantified as expected.
% Must run above func_detect_measure_bacteria section prior to this.

% Specify the condition folder and the individual image you want to test (file name without '_BF.png'):
test_condition = 'Full Template';
test_sample = 'B.1. Full Template';
% restrict range of measured fluorescence values:
minyvalue = 0; % define minimum value
maxyvalue = 2 * positive_threshold; % define maximum value, default is twice the GFP positive gate CHANGE THIS

% plot the image
f = find(strcmp(condition_folders, test_condition), 1);
func_plot_MES_for_sample_bacteria(test_condition, test_sample, allMES{f}, image_names{f}, adapt_sensitivity, min_cell_area, max_cell_area, maxyvalue, minyvalue)

%% ----------------- MAKE % FLUORESCENT CELLS VIOLIN PLOT, IMAGES -------------------------------------------------------------------------------------
% plot the percentage of GFP positive cells for each condition
%       - each dot is one image (the average of the images in that repeat)
%       - boxplots indicate the interquartile range
%       - notches indicate the median
%       - asterisks indicate the mean

% % GFP positive cells for each repeat
percentPositive = cell(1, length(condition_folders));
for f = 1:length(condition_folders)
    percentPositive{f} = repeatTable.PercentPositive(strcmp(repeatTable.Condition, condition_folders{f}));
end

% calculate % GFP positive cells for each image
percentPositive = cell(1, length(condition_folders));
for f = 1:length(condition_folders)
    percentPositive{f} = percentTable.PercentPositive(strcmp(percentTable.Condition, condition_folders{f}));
end

% Create violin plot grouped by condition
figure3 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);

for j = 1:length(percentPositive)
    if ~isempty(percentPositive{j})
        func_al_goodplot(percentPositive{j}, j*n, [], manualColors{j});
    end
end
%set(findobj(gca, 'Type', 'scatter'), 'SizeData', 60) % bigger dots, as there is only one per repeat

ylabel('GFP-expressing bacteria (% of total)') % define y axis label
ylim([0 60]) % COMMENT for automatic defined y axis limits

xlim([0 (length(percentPositive)*n)+1])
xticks(n * (1:length(percentPositive)))
%xticklabels(cleanNames);
xticklabels({'Inactive Template', 'Full Template', 'MagSNA(-AMF)', 'MagSNA(+AMF)'});
xtickangle(30)
grid off;

ax1 = gca;
ax1.YAxis.FontSize = 18;
ax1.XAxis.FontSize = 18;
ax1.YAxis.Exponent = 0;
ax1.LineWidth=2;
box off

% NOW RUN THIS SECTION ^

%% ----------------- MAKE % FLUORESCENT CELLS VIOLIN PLOT, IMAGES, FOLD-CHANGE -------------------------------------------------------------------------------------
% plot the percentage of GFP positive cells for each condition
%       - each dot is one image (the average of the images in that repeat)
%       - boxplots indicate the interquartile range
%       - notches indicate the median
%       - asterisks indicate the mean

% calculate % GFP positive cells for each image
percentPositive = cell(1, length(condition_folders));
for f = 1:length(condition_folders)
    percentPositive{f} = percentTable.PercentPositive(strcmp(percentTable.Condition, condition_folders{f}));
end

% convert to fold change relative to the negative control
negMean = mean(percentPositive{isNegative});
foldChange = cellfun(@(x) x / negMean, percentPositive, 'UniformOutput', false);

% Create violin plot grouped by condition
figure3 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);

for j = 1:length(foldChange)
    if ~isempty(foldChange{j})
        func_al_goodplot(foldChange{j}, j*n, [], manualColors{j});
    end
end
%set(findobj(gca, 'Type', 'scatter'), 'SizeData', 60) % bigger dots, as there is only one per repeat

ylabel({'GFP-expressing cells', '(fold change vs Inactive Template)'})
ylim([0 8]) % COMMENT for automatic defined y axis limits

xlim([0 (length(foldChange)*n)+1])
xticks(n * (1:length(foldChange)))
%xticklabels(cleanNames);
xticklabels({'Inactive Template', 'Full Template', 'MagSNA(-AMF)', 'MagSNA(+AMF)'});
xtickangle(30)
grid off;

ax1 = gca;
ax1.YAxis.FontSize = 18;
ax1.XAxis.FontSize = 18;
ax1.YAxis.Exponent = 0;
ax1.LineWidth=2;
box off
% NOW RUN THIS SECTION ^

%% ----------------- SWITCH RESPONSE, NORMALISED TO THE CONTROLS, IMAGES -------------------------------------------------------------------------------------
% 0% = average of the negative control images (Inactive Template), 100% = average of the positive control images (Full Template)
% each dot = one image, black line & error bar = mean ± SD

show = condition_folders(ismember(condition_folders, {'-AMF', '+AMF'})); % conditions to plot, use show = condition_folders; to include the controls

pct = @(name) percentTable.PercentPositive(strcmp(percentTable.Condition, name) & ~isnan(percentTable.PercentPositive));
negMean = mean(pct('Inactive Template'));
posMean = mean(pct('Full Template'));

figure6 = figure('Renderer', 'painters', 'Position', [10 10 400 600]);
hold on
for k = 1:numel(show)
    values = (pct(show{k}) - negMean) / (posMean - negMean) * 100;
    scatter(k + 0.2 * (rand(size(values)) - 0.5), values, 60, manualColors{strcmp(condition_folders, show{k})}, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7)
    plot(k + [-0.25 0.25], mean(values) * [1 1], 'k', 'LineWidth', 2)
    errorbar(k, mean(values), std(values), 'k', 'LineWidth', 1.5, 'CapSize', 12)
end
%yline(0, ':k', 'Inactive Template', 'FontSize', 12); yline(100, ':k', 'Full Template', 'FontSize', 12);
yline(0, '--'); yline(100, '--');
xlim([0.5 numel(show) + 0.5]); xticks(1:numel(show)); xticklabels(strrep(show, '_', '\_'));
ylabel({'GFP positive cells', '(% of Full Template response)'});
ylim([-100 300])
ax = gca; ax.FontSize = 16; ax.LineWidth = 2; box off

%% ----------------- MAKE % FLUORESCENT CELLS VIOLIN PLOT, REPEATS -------------------------------------------------------------------------------------
% plot the percentage of GFP positive cells for each condition
%       - each dot is one repeat (the average of the images in that repeat)
%       - boxplots indicate the interquartile range
%       - notches indicate the median
%       - asterisks indicate the mean

% % GFP positive cells for each repeat
percentPositive = cell(1, length(condition_folders));
for f = 1:length(condition_folders)
    percentPositive{f} = repeatTable.PercentPositive(strcmp(repeatTable.Condition, condition_folders{f}));
end

% Create violin plot grouped by condition
figure3 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);

for j = 1:length(percentPositive)
    if ~isempty(percentPositive{j})
        func_al_goodplot(percentPositive{j}, j*n, [], manualColors{j});
    end
end
set(findobj(gca, 'Type', 'scatter'), 'SizeData', 60) % bigger dots, as there is only one per repeat

ylabel('GFP-expressing bacteria (% of total)') % define y axis label
ylim([0 60]) % COMMENT for automatic defined y axis limits

xlim([0 (length(percentPositive)*n)+1])
xticks(n * (1:length(percentPositive)))
%xticklabels(cleanNames);
xticklabels({'Inactive Template', 'Full Template', 'MagSNA(-AMF)', 'MagSNA(+AMF)'});
xtickangle(30)
grid off;

ax1 = gca;
ax1.YAxis.FontSize = 18;
ax1.XAxis.FontSize = 18;
ax1.YAxis.Exponent = 0;
ax1.LineWidth=2;
box off

% NOW RUN THIS SECTION ^

%% ----------------- MAKE % FLUORESCENT CELLS VIOLIN PLOT, REPEATS, FOLD CHANGE -------------------------------------------------------------------------------------
% plot the percentage of GFP positive cells, relative to the negative control (mean of the Inactive Template repeats = 1)
%       - each dot is one repeat (the average of the images in that repeat)
%       - boxplots indicate the interquartile range
%       - notches indicate the median
%       - asterisks indicate the mean

% % GFP positive cells for each repeat
percentPositive = cell(1, length(condition_folders));
for f = 1:length(condition_folders)
    percentPositive{f} = repeatTable.PercentPositive(strcmp(repeatTable.Condition, condition_folders{f}));
end

% fold change relative to each repeat's own negative control
foldChange = cell(1, length(condition_folders));
for f = 1:length(condition_folders)
    values = [];
    for r = 1:numel(repeats)
        neg = repeatTable.PercentPositive(strcmp(repeatTable.Condition, negative_control) & strcmp(repeatTable.Repeat, repeats{r}));
        x = repeatTable.PercentPositive(strcmp(repeatTable.Condition, condition_folders{f}) & strcmp(repeatTable.Repeat, repeats{r}));
        if isempty(neg) || isempty(x), continue, end % skip repeats missing the negative control
        values(end+1) = x / neg; %#ok<SAGROW>
    end
    foldChange{f} = values;
end

% Create violin plot grouped by condition
figure5 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);

for j = 1:length(foldChange)
    if ~isempty(foldChange{j})
        func_al_goodplot(foldChange{j}, j*n, [], manualColors{j});
    end
end
set(findobj(gca, 'Type', 'scatter'), 'SizeData', 60) % bigger dots, as there is only one per repeat

ylabel({'GFP-expressing cells', '(fold change vs Inactive Template)'})
ylim([0 7]) % COMMENT for automatic defined y axis limits

xlim([0 (length(foldChange)*n)+1])
xticks(n * (1:length(foldChange)))
%xticklabels(cleanNames);
xticklabels({'Inactive Template', 'Full Template', 'MagSNA(-AMF)', 'MagSNA(+AMF)'});
xtickangle(30)
grid off;

ax1 = gca;
ax1.YAxis.FontSize = 18;
ax1.XAxis.FontSize = 18;
ax1.YAxis.Exponent = 0;
ax1.LineWidth=2;
box off

% NOW RUN THIS SECTION ^

%% ----------------- SUPERPLOT OF % GFP POSITIVE CELLS -------------------------------------------------------------------------------------
% small dots = individual images, coloured by repeat
% large dots = mean of each repeat
% black line & error bar = mean ± SD of the repeat means

repeatColors = [0.00 0.45 0.70;   % repeat A, blue
                0.84 0.37 0.00;   % repeat B, dark orange
                0.00 0.62 0.45];  % repeat C, green

figure4 = figure('Renderer', 'painters', 'Position', [10 10 500 700]);
hold on
for c = 1:length(condition_folders)
    inCondition = strcmp(percentTable.Condition, condition_folders{c}) & ~isnan(percentTable.PercentPositive);
    repeatMeans = [];
    for r = 1:numel(repeats)
        rows = inCondition & strcmp(percentTable.Repeat, repeats{r});
        if ~any(rows), continue, end
        xRepeat = c + 0.4 * ((r - 0.5) / numel(repeats) - 0.5); % spread repeats side by side
        scatter(xRepeat + 0.04 * randn(sum(rows), 1), percentTable.PercentPositive(rows), 30, repeatColors(r,:), 'filled', 'MarkerFaceAlpha', 0.4)
        repeatMeans(end+1) = mean(percentTable.PercentPositive(rows)); %#ok<SAGROW>
        scatter(xRepeat, repeatMeans(end), 50, repeatColors(r,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5)
    end
    plot(c + [-0.3 0.3], mean(repeatMeans) * [1 1], 'k', 'LineWidth', 3)
    errorbar(c, mean(repeatMeans), std(repeatMeans), 'k', 'LineWidth', 1.5, 'CapSize', 12)
end

ylabel('GFP Positive Cells (%)')
xlim([0.5 length(condition_folders) + 0.5])
xticks(1:length(condition_folders)); xticklabels(cleanNames); xtickangle(30)
ax1 = gca; ax1.FontSize = 18; ax1.LineWidth = 2; box off

%% ----------------- SAVE FIGURES, TABLES & PARAMETERS -------------------------------------------------------------------------------------
% Saves every open figure as a 300 dpi png, the % positive tables as spreadsheets,
% and the parameters used, into a new time-stamped folder so results are never overwritten.

output_folder = fullfile(working_folder, 'results', char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss')));
mkdir(output_folder);

figs = findobj('Type', 'figure');
for k = 1:numel(figs)
    exportgraphics(figs(k), fullfile(output_folder, sprintf('figure_%d.png', figs(k).Number)), 'Resolution', 300);
end

writetable(percentTable, fullfile(output_folder, 'percent_positive_per_image.csv'));
writetable(repeatTable, fullfile(output_folder, 'percent_positive_per_repeat.csv'));

parameters = table(adapt_sensitivity, min_cell_area, max_cell_area, string(negative_control), gate_percent, positive_threshold, scale_factor, ...
    'VariableNames', {'adapt_sensitivity', 'min_cell_area', 'max_cell_area', 'negative_control', 'gate_percent', 'positive_threshold', 'scale_factor'});
writetable(parameters, fullfile(output_folder, 'parameters.csv'));

fprintf('Saved %d figures, the %% positive tables and the parameters to %s\n', numel(figs), output_folder);
