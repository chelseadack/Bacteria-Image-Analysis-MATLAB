function [image_names, image_repeats, grey_DYE_images, grey_MES_images] = func_load_images_bacteria(condition_folder)
% FUNC_LOAD_IMAGES_BACTERIA  Load all of the image files for one condition folder.
%
%   [image_names, image_repeats, grey_DYE_images, grey_MES_images] = func_load_images_bacteria(condition_folder)
%
%   Data must be in this format:
%       - One folder per condition, containing all of the repeats for that condition.
%       - Colour or greyscale images must be saved as '.png'.
%       - Each unique image must have 2 images named 'repeat.image. condition':
%           Eg. 'B.3. -AMF_BF.png'  = Brightfield image used to detect the cells (repeat B, image 3).
%               'B.3. -AMF_GFP.png' = Image of fluorescence to measure.
%
%   Input parameters:
%       - condition_folder = Name of the condition folder eg. '-AMF'.
%   Output (all 1xN cell arrays, one entry per image):
%       - image_names = Image names, without '_BF.png' eg. 'B.3. -AMF'.
%       - image_repeats = Repeat each image belongs to (the text before the first '.') eg. 'B'.
%       - grey_DYE_images = Greyscale brightfield images.
%       - grey_MES_images = Greyscale GFP images.
%
%   See also FUNC_DETECT_MEASURE_BACTERIA.

    arguments
        condition_folder {mustBeTextScalar}
    end
    if ~isfolder(condition_folder)
        error('func_load_images_bacteria:noFolder', 'Condition folder "%s" not found in %s', condition_folder, pwd);
    end

    fprintf('Started loading images from %s\n', condition_folder);

    files = dir(fullfile(condition_folder, '*_BF.png')); % Find all brightfield images in the folder
    nfiles = length(files); % Number of images found
    if nfiles == 0
        warning('func_load_images_bacteria:noImages', 'No *_BF.png images found in %s', condition_folder);
    end

    image_names = cell(1, nfiles);
    image_repeats = cell(1, nfiles);
    grey_DYE_images = cell(1, nfiles);
    grey_MES_images = cell(1, nfiles);

    % For each brightfield image, load it and its matching GFP image
    for i = 1:nfiles

        image_names{i} = erase(files(i).name, '_BF.png');

        % the repeat is the text before the first '.', eg. 'B' from 'B.3. -AMF'
        tokens = regexp(image_names{i}, '^([^.]+)\.[^.]+\.', 'tokens', 'once');
        if isempty(tokens)
            error('func_load_images_bacteria:badName', 'Image "%s" does not follow the naming pattern ''repeat.image. condition_BF.png''', fullfile(condition_folder, files(i).name));
        end
        image_repeats{i} = tokens{1};

        MES_filename = fullfile(condition_folder, [image_names{i} '_GFP.png']);
        if ~isfile(MES_filename)
            error('func_load_images_bacteria:noGFP', 'No matching GFP image found for %s', fullfile(condition_folder, files(i).name));
        end

        % Converts the files to greyscale and groups data into respective array
        grey_DYE_images{i} = im2gray(imread(fullfile(condition_folder, files(i).name)));
        grey_MES_images{i} = im2gray(imread(MES_filename));

    end

    fprintf('Finished loading %d images from %s\n', nfiles, condition_folder);

end
