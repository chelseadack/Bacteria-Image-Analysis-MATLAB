function func_test_detection_bacteria(test_condition, test_sample, adapt_sensitivity, min_cell_area, max_cell_area)
% FUNC_TEST_DETECTION_BACTERIA  Show where cells are detected in one image, to check the detection parameters.
%
%   func_test_detection_bacteria(test_condition, test_sample, adapt_sensitivity, min_cell_area, max_cell_area)
%
%   For the specified test sample, it loads the images then creates a figure showing where it detected cells
%   and where it will measure fluorescence.
%
%   Input parameters:
%       - test_condition = Name of the condition folder containing the test sample eg. '-AMF'.
%       - test_sample = Name of specific sample to test eg. '1.1 -AMF'.
%           - Test sample images must have these naming patterns:
%               Eg. '1.1 -AMF_BF.png'  = Brightfield image used to detect the cells.
%                   '1.1 -AMF_GFP.png' = Image of fluorescence to measure.
%       - adapt_sensitivity = Noise factor for the threshold (higher = fewer pixels called cell), see func_detect_bacteria.
%       - min_cell_area = Minimum cell area to detect in pixels.
%       - max_cell_area = Maximum cell area to detect in pixels.
%   Output:
%       - Brightfield image with overlay of where it detected cells in the image.
%       - GFP image with overlay of where it will measure fluorescence in the image.
%
%   See also FUNC_DETECT_BACTERIA.

    fprintf('Started testing cell detection parameters\n');

    %%%%%%%%%%%%% LOAD IMAGES FOR SPECIFIED SAMPLE %%%%%%%%%%%%

    % load images for the individual sample you specified above
    BF_filename = fullfile(test_condition, [test_sample '_BF.png']);
    GFP_filename = fullfile(test_condition, [test_sample '_GFP.png']);
    if ~isfile(BF_filename) || ~isfile(GFP_filename)
        error('func_test_detection_bacteria:noImage', 'Could not find %s and %s in %s', BF_filename, GFP_filename, pwd);
    end
    markerImage = imread(BF_filename);
    detectionImage = imread(GFP_filename);

    % convert the brightfield image to greyscale
    greyMarkerImage = im2gray(markerImage);

    % Find the cells in the brightfield image
    [cells, cellMask] = func_detect_bacteria(greyMarkerImage, adapt_sensitivity, min_cell_area, max_cell_area);

    %%%%%%%%%%%%%%% PLOT FIGURE OF DETECTED CELLS %%%%%%%%%%%%%

    % this will create a figure showing where it detected cells and where it will measure fluorescence.

    figure('Renderer', 'painters', 'Position', [50 100 1400 600]);
    subplot(1,2,1)

    imshow(adapthisteq(greyMarkerImage)) % contrast enhanced so cells are easier to see by eye
    hold on
    visboundaries(cellMask, 'Color', 'c', 'LineWidth', 0.2);
    text(50,-100,sprintf('%d cells identified in brightfield image', numel(cells)),'Fontsize',18)
    hold off

    % shows the location of the cells on the fluorescence image
    subplot(1,2,2)
    imshow(detectionImage, []) % [] = contrast stretched, as GFP images are very dim
    hold on % keep the underlying image visible
    visboundaries(cellMask, 'Color', 'c', 'LineWidth', 0.2);
    text(80,-100,'GFP expression in those cells','Fontsize',18)
    hold off

    fprintf('Finished testing cell detection parameters\n');

end
