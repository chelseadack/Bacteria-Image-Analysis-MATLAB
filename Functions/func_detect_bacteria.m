function [cells, cellMask] = func_detect_bacteria(grey_DYE_image, adapt_sensitivity, min_cell_area, max_cell_area)
% FUNC_DETECT_BACTERIA  Find the outline of every bacterial cell in a brightfield image.
%
%   [cells, cellMask] = func_detect_bacteria(grey_DYE_image, adapt_sensitivity, min_cell_area, max_cell_area)
%
%   This function is shared by func_test_detection_bacteria, func_detect_measure_bacteria and
%   func_plot_MES_for_sample_bacteria, so all three always detect exactly the same cells, in the same order.
%
%   Method:
%       1. Flatten uneven illumination by subtracting a heavily blurred copy of the image.
%       2. Cells can be darker OR brighter than the background, so use the size of the difference (contrast image).
%       3. Threshold the contrast image at two levels, both relative to a robust (median / MAD) noise estimate:
%           - high threshold: noise estimated from the whole image. Works best inside clusters, leaves gaps between cells.
%           - low threshold: noise estimated from background pixels only. Works best for free cells.
%          Use the low threshold for free cells and the high threshold inside clusters.
%       4. Split touching cells with a distance-transform watershed.
%       5. Keep only objects within the defined size range.
%
%   Input parameters:
%       - grey_DYE_image = Greyscale brightfield image used to detect the cells.
%       - adapt_sensitivity = Noise factor for the threshold: a pixel is called cell if it is this many noise levels
%                             above the background (higher = fewer pixels called cell, lower = more pixels called cell).
%       - min_cell_area = Minimum cell area to keep, in pixels.
%       - max_cell_area = Maximum cell area to keep, in pixels.
%   Output:
%       - cells = structure array with one entry per detected cell ('Area', 'Centroid', 'PixelIdxList').
%       - cellMask = logical image of all detected cells.
%
%   Requires the Image Processing Toolbox.
%
%   See also FUNC_DETECT_MEASURE_BACTERIA, FUNC_TEST_DETECTION_BACTERIA.

    arguments
        grey_DYE_image (:,:) {mustBeNumeric}
        adapt_sensitivity (1,1) {mustBeNumeric, mustBePositive}
        min_cell_area (1,1) {mustBeNumeric, mustBeNonnegative}
        max_cell_area (1,1) {mustBeNumeric, mustBePositive}
    end
    if min_cell_area > max_cell_area
        error('func_detect_bacteria:areaRange', 'min_cell_area (%g) must not be larger than max_cell_area (%g).', min_cell_area, max_cell_area);
    end

    % Fixed detection constants
    illumination_sigma = 50; % blur (pixels) used to estimate uneven illumination
    noise_sigma = 1; % slight blur (pixels) to reduce speckle noise
    mad_to_sd = 1.4826; % converts median absolute deviation to standard deviation for normally distributed noise
    cluster_area = 3 * max_cell_area; % objects bigger than this (pixels) in the low threshold mask are treated as clusters
    split_depth = 2; % ignore dips in the distance transform shallower than this (pixels) when splitting touching cells

    img = im2double(grey_DYE_image); % grey_DYE_image is uint8 which can't be < 0, so convert to detect dark cells too

    % Flatten uneven illumination
    flatImage = img - imgaussfilt(img, illumination_sigma);

    % Cells can be darker or brighter than background, so use the size of the difference
    contrastImage = imgaussfilt(abs(flatImage), noise_sigma); % slight blur reduces speckle noise

    % High threshold: noise estimated from the whole image (works best inside clusters, leaves gaps between cells)
    backgroundLevel = median(contrastImage(:));
    noiseLevel = mad_to_sd * median(abs(contrastImage(:) - backgroundLevel)); % robust estimate of noise
    highMask = contrastImage > backgroundLevel + adapt_sensitivity * noiseLevel;
    highMask = imopen(highMask, strel('disk', 1)); % removes leftover single-pixel speckle

    % Low threshold: noise estimated from background pixels only (works best for free cells)
    backgroundPixels = contrastImage(:);
    for iteration = 1:5
        backgroundLevel = median(backgroundPixels);
        noiseLevel = mad_to_sd * median(abs(backgroundPixels - backgroundLevel));
        backgroundPixels = backgroundPixels(backgroundPixels < backgroundLevel + 3 * noiseLevel); % remove cell pixels & re-estimate
    end
    lowMask = contrastImage > backgroundLevel + adapt_sensitivity * noiseLevel;
    lowMask = imopen(lowMask, strel('disk', 1)); % removes leftover single-pixel speckle
    lowMask = imfill(lowMask, 'holes'); % fills pale cell centres inside dark outlines

    % Use the low threshold for free cells, but the high threshold inside clusters
    clusters = bwareafilt(lowMask, [cluster_area Inf], 4);
    cellMask = (lowMask & ~clusters) | (highMask & clusters);

    % Distance transform and watershed to cut lines between touching cells
    D = -bwdist(~cellMask);
    D = imhmin(D, split_depth); % ignore dips shallower than split_depth pixels
    D(~cellMask) = -Inf;
    L = watershed(D);
    cellMask(L == 0) = 0;

    % Label each cell (4-connectivity so cells split by the watershed stay separate)
    cells = regionprops(bwconncomp(cellMask, 4), 'Area', 'Centroid', 'PixelIdxList');

    % Keep only cells within the defined size range
    if ~isempty(cells)
        cells = cells([cells.Area] >= min_cell_area & [cells.Area] <= max_cell_area);
    end

    % Remake the mask so it only contains the kept cells
    cellMask = false(size(cellMask));
    cellMask(vertcat(cells.PixelIdxList)) = true;

end
