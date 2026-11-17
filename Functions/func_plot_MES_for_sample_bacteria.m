function func_plot_MES_for_sample_bacteria(test_condition, test_sample, allMES, image_names, adapt_sensitivity, min_cell_area, max_cell_area, maxyvalue, minyvalue)
% FUNC_PLOT_MES_FOR_SAMPLE_BACTERIA  Overlay each cell's measured fluorescence onto its GFP image.
%
%   func_plot_MES_for_sample_bacteria(test_condition, test_sample, allMES, image_names, ...
%                                     adapt_sensitivity, min_cell_area, max_cell_area, maxyvalue, minyvalue)
%
%   Used to check the cells are being quantified as expected. Must run func_detect_measure_bacteria first.
%
%   Input parameters:
%       - test_condition = name of the condition folder containing the sample (e.g., '-AMF').
%       - test_sample = name of the sample (e.g., '1.1 -AMF').
%       - allMES = cell array returned by func_detect_measure_bacteria for this condition.
%       - image_names = image names returned by func_load_images_bacteria for this condition.
%       - adapt_sensitivity = noise factor for the threshold used in detection.
%       - min_cell_area = minimum cell area in pixels.
%       - max_cell_area = maximum cell area in pixels.
%       - maxyvalue = maximum measured fluorescence value to plot.
%       - minyvalue = minimum measured fluorescence value to plot.
%   Output:
%       - Overlays a colour circle indicating measured fluorescence value onto each cell in the test image.
%
%   See also FUNC_DETECT_BACTERIA, FUNC_DETECT_MEASURE_BACTERIA.

    if minyvalue >= maxyvalue
        error('func_plot_MES_for_sample_bacteria:range', 'minyvalue (%g) must be smaller than maxyvalue (%g).', minyvalue, maxyvalue);
    end

    %%%%%%%%%%%%% LOAD IMAGES FOR SPECIFIED SAMPLE & DETECT CELLS %%%%%%%%%%%%%

    % load images for the individual sample you specified above
    markerImage = imread(fullfile(test_condition, [test_sample '_BF.png']));
    detectionImage = imread(fullfile(test_condition, [test_sample '_GFP.png']));

    % convert the brightfield image to greyscale
    greyMarkerImage = im2gray(markerImage);

    % Find the cells in the brightfield image (same detection, so same order as allMES)
    cells = func_detect_bacteria(greyMarkerImage, adapt_sensitivity, min_cell_area, max_cell_area);
    centers = vertcat(cells.Centroid);

    %%%%%%%%%%%%%%% FIND ASSIGNED VALUES FOR THIS SAMPLE %%%%%%%%%%%%%%%%%

    % Find the index of this sample within the condition
    j = find(strcmp(image_names, test_sample), 1);
    if isempty(j)
        error('func_plot_MES_for_sample_bacteria:noSample', 'Sample %s not found in condition %s', test_sample, test_condition);
    end

    MES_values = allMES{j};

    % The cells are matched to their values by order, so the detection must be the same as when allMES was made
    if numel(MES_values) ~= numel(cells)
        error('func_plot_MES_for_sample_bacteria:mismatch', ...
            '%d cells detected but %d values stored for %s. Re-run the detect & measure section with the current parameters.', ...
            numel(cells), numel(MES_values), test_sample);
    end

    %%%%%%%%%%%%% OVERLAY COLOUR-CODED VALUES ONTO IMAGE %%%%%%%%%%%%%

    figure('Renderer','painters','Position',[50 50 1400 600],'Visible','on');

    subplot(1,2,1)
    imshow(detectionImage, []) % [] = contrast stretched, as GFP images are very dim

    subplot(1,2,2)
    imshow(detectionImage, []) % same image, with the measured values overlaid
    hold on

    % Make colour-coded MES spots
    cmap = turbo(256); % colour scheme
    vmin = minyvalue; % minimum value
    vmax = maxyvalue; % maximum value
    MES_clipped = min(max(MES_values, vmin), vmax);
    MES_scaled = round( (MES_clipped - vmin) / (vmax - vmin) * 255 ) + 1;
    MES_scaled(MES_scaled<1) = 1;
    MES_scaled(MES_scaled>256) = 256;
    spot_colors = cmap(MES_scaled,:);

    % plot coloured spots onto image
    for k = 1:numel(MES_values)
        plot(centers(k,1), centers(k,2), 'o', 'MarkerFaceColor', spot_colors(k,:), ...
             'MarkerEdgeColor','k','MarkerSize',6);
    end

    title(['Per-cell MES values: ', test_sample])
    hold off

    % Separate colorbar figure
    figure('Position',[1200 100 100 400]);
    ax = axes;
    colormap(ax, turbo(256))
    clim(ax, [vmin vmax])
    cb = colorbar(ax);
    ylabel(cb,'MES value')
    cb.Position(1) = 0.38;
    cb.Position(3) = 0.25;
    axis off

end
