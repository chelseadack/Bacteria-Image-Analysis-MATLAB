function allMES = func_detect_measure_bacteria(adapt_sensitivity, min_cell_area, max_cell_area, grey_DYE_images, grey_MES_images)
% FUNC_DETECT_MEASURE_BACTERIA  Detect the bacteria in each image, then measure average whole cell fluorescence.
%
%   allMES = func_detect_measure_bacteria(adapt_sensitivity, min_cell_area, max_cell_area, grey_DYE_images, grey_MES_images)
%
%   No background is subtracted: GFP positive cells are found by gating against the negative control instead.
%
%   Input parameters:
%       - adapt_sensitivity = Noise factor for the threshold (higher = fewer pixels called cell), see func_detect_bacteria.
%       - min_cell_area = Minimum cell area to detect in pixels.
%       - max_cell_area = Maximum cell area to detect in pixels.
%       - grey_DYE_images = Collection of all brightfield images loaded using func_load_images_bacteria().
%       - grey_MES_images = Collection of all GFP images loaded using func_load_images_bacteria().
%   Output:
%       - allMES = 1xN cell array, one entry per image. Each entry is a row vector of the mean GFP fluorescence
%                  of every cell detected in that image (same order as func_detect_bacteria).
%
%   See also FUNC_DETECT_BACTERIA, FUNC_LOAD_IMAGES_BACTERIA.

    arguments
        adapt_sensitivity (1,1) {mustBeNumeric, mustBePositive}
        min_cell_area (1,1) {mustBeNumeric, mustBeNonnegative}
        max_cell_area (1,1) {mustBeNumeric, mustBePositive}
        grey_DYE_images cell
        grey_MES_images cell
    end
    if numel(grey_DYE_images) ~= numel(grey_MES_images)
        error('func_detect_measure_bacteria:sizeMismatch', 'There must be the same number of brightfield and GFP images.');
    end

    fprintf('Started cell detection\n');

    % create empty arrays
    allMES = cell(1, length(grey_DYE_images));

    % for each sample, detect cells and measure average fluorescence
    for j = 1:length(grey_DYE_images)

        % find cells
        cells = func_detect_bacteria(grey_DYE_images{j}, adapt_sensitivity, min_cell_area, max_cell_area);

        img = grey_MES_images{j};

        %%%%%%%%%%%%%% MEASURE WHOLE CELL %%%%%%%%%%%%%

        % for every cell detected, calculate average whole cell fluorescence
        allMES{j} = zeros(1, numel(cells));
        for k = 1:numel(cells)
            allMES{j}(k) = mean(img(cells(k).PixelIdxList));
        end

        fprintf('    image %d of %d: %d cells\n', j, length(grey_DYE_images), numel(cells));

    end

    fprintf('Finished cell detection\n');

end
