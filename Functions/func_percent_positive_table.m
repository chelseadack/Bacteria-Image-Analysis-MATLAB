function [percentTable, repeatTable] = func_percent_positive_table(condition_folders, image_names, image_repeats, allMES, repeats, positive_threshold)
% FUNC_PERCENT_POSITIVE_TABLE  % GFP positive cells in every image, and the average for every repeat.
%
%   [percentTable, repeatTable] = func_percent_positive_table(condition_folders, image_names, image_repeats, allMES, repeats, positive_threshold)
%
%   Every image is gated with the same threshold.
%
%   Input parameters:
%       - condition_folders = 1xC cell array of condition names.
%       - image_names = 1xC cell array, image names for each condition (from func_load_images_bacteria).
%       - image_repeats = 1xC cell array, repeat of each image for each condition (from func_load_images_bacteria).
%       - allMES = 1xC cell array, per-image fluorescence values for each condition (from func_detect_measure_bacteria).
%       - repeats = 1xR cell array of repeat names eg. {'A', 'B', 'C'}.
%       - positive_threshold = GFP positive gate (from func_gate_threshold).
%   Output:
%       - percentTable = one row per image, columns Condition, Repeat, Image, TotalCells, PositiveCells, PercentPositive.
%                        PercentPositive is NaN for images with no cells detected.
%       - repeatTable = one row per repeat of each condition, columns Condition, Repeat, Images, TotalCells, PercentPositive,
%                       where PercentPositive is the average of that repeat's images (images with no cells are left out).

    % one row per image
    Condition = {}; Repeat = {}; Image = {}; TotalCells = []; PositiveCells = [];
    for f = 1:length(condition_folders)
        for i = 1:length(image_names{f})
            Condition{end+1,1} = condition_folders{f}; %#ok<AGROW>
            Repeat{end+1,1} = image_repeats{f}{i}; %#ok<AGROW>
            Image{end+1,1} = image_names{f}{i}; %#ok<AGROW>
            TotalCells(end+1,1) = numel(allMES{f}{i}); %#ok<AGROW>
            PositiveCells(end+1,1) = sum(allMES{f}{i} > positive_threshold); %#ok<AGROW>
        end
    end
    PercentPositive = PositiveCells ./ TotalCells * 100;
    percentTable = table(Condition, Repeat, Image, TotalCells, PositiveCells, PercentPositive);

    % one row per repeat of each condition, average of its images
    Condition = {}; Repeat = {}; Images = []; TotalCells = []; PercentPositive = [];
    for f = 1:length(condition_folders)
        for r = 1:length(repeats)
            rows = strcmp(percentTable.Condition, condition_folders{f}) & strcmp(percentTable.Repeat, repeats{r}) & ~isnan(percentTable.PercentPositive);
            if ~any(rows), continue, end
            Condition{end+1,1} = condition_folders{f}; %#ok<AGROW>
            Repeat{end+1,1} = repeats{r}; %#ok<AGROW>
            Images(end+1,1) = sum(rows); %#ok<AGROW>
            TotalCells(end+1,1) = sum(percentTable.TotalCells(rows)); %#ok<AGROW>
            PercentPositive(end+1,1) = mean(percentTable.PercentPositive(rows)); %#ok<AGROW>
        end
    end
    repeatTable = table(Condition, Repeat, Images, TotalCells, PercentPositive);

end
