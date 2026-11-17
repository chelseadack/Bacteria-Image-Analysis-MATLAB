# Bacteria GFP quantification (MATLAB)

Automated detection of bacteria in 100× brightfield microscope images, with per-cell GFP fluorescence measurement and population-level analysis of a genetic switch.

The pipeline finds every cell in a brightfield image, measures the mean GFP fluorescence in each one, and compares conditions using:

- per-cell distributions on an asinh scale (swarm plots and histograms), for each biological repeat and for the whole population
- a single GFP-positive gate set from the negative control cells of all repeats (95th percentile), used for every repeat. No background subtraction is needed, because cells are classified relative to the negative control.
- the % of GFP-positive cells and the fold change vs the negative control, with one data point per biological repeat
- a switch response normalised between the negative (0%) and positive (100%) controls

## Cell detection

Bacteria in brightfield can appear darker **or** brighter than the background, and they are often packed into dense aggregates. `func_detect_bacteria.m` handles both:

1. **Flatten illumination:** subtract a Gaussian-blurred (σ = 50 px) copy of the image.
2. **Contrast image:** take the absolute deviation from the background, so both dark and bright cells become bright objects.
3. **Dual robust threshold:** the noise is estimated with the median absolute deviation (MAD), which is not biased by the cells themselves.
   - The *low threshold* estimates noise from background pixels only, with iterative outlier rejection. It segments free cells well.
   - The *high threshold* estimates noise from the whole image. It separates cells inside dense clusters.
   - Objects larger than 3 × the maximum cell area are treated as clusters and use the high-threshold mask. Everything else uses the low-threshold mask.
4. **Split touching cells** with a distance-transform watershed (h-minima suppression).
5. **Size filter:** keep objects between `min_cell_area` and `max_cell_area`.

The same function is called by every step that needs cell outlines, so the cells are always detected identically and in the same order.

## Repository structure

```Bacteria-Image-Analysis-MATLAB/
├── MAIN_SCRIPT_BACTERIA.m      # main analysis script, run it section by section
├── Functions/                  # functions called by the main script
└── <condition folders>         # your images (see Data format)
```

| File | Purpose |
|---|---|
| `MAIN_SCRIPT_BACTERIA.m` | Main analysis script. Run it section by section. |
| `Functions/func_load_images_bacteria.m` | Loads the brightfield/GFP image pairs for one condition folder |
| `Functions/func_detect_bacteria.m` | Cell segmentation (see above) |
| `Functions/func_detect_measure_bacteria.m` | Detects cells in every image and measures the mean GFP per cell |
| `Functions/func_gate_threshold.m` | GFP-positive gate (percentile of the negative control cells) |
| `Functions/func_percent_positive_table.m` | % GFP-positive cells per image, and the average of each repeat |
| `Functions/func_plot_swarm_bacteria.m` | Swarm plot of every cell, with the gate |
| `Functions/func_plot_smoothed_histogram_bacteria.m` | Smoothed histograms of every cell, with the gate |
| `Functions/func_test_detection_bacteria.m` | QC figure: detected cell outlines on the brightfield and GFP images |
| `Functions/func_plot_MES_for_sample_bacteria.m` | QC figure: each cell coloured by its measured fluorescence |
| `Functions/func_al_goodplot.m` | Third-party violin/box plot function (see Acknowledgements) |

## Data format

```
working_folder/
├── MAIN_SCRIPT_BACTERIA.m
├── Functions/
├── Inactive Template/          # one folder per condition
│   ├── B.1. Inactive Template_BF.png   # repeat B, image 1, brightfield (used for detection)
│   ├── B.1. Inactive Template_GFP.png  # repeat B, image 1, fluorescence (measured)
│   └── ...
├── Full Template/
├── -AMF/
└── +AMF/
```

- The images can be colour or greyscale `.png` files.
- Images are named `repeat.image. condition`. The text before the first `.` is the biological repeat.

## How to run

1. Requirements: MATLAB R2022a or later with the **Image Processing Toolbox**. It has been tested in MATLAB Online.
2. Put `MAIN_SCRIPT_BACTERIA.m`, the `Functions/` folder and the condition folders in one folder (`working_folder`). The script adds `Functions/` to the MATLAB path.
3. Open `MAIN_SCRIPT_BACTERIA.m` and set `working_folder`, `condition_folders` and the gate settings (`negative_control`, `gate_percent`).
4. Run the sections in order (Ctrl/Cmd + Enter):
   1. **Input parameters:** the detection settings (`adapt_sensitivity`, `min_cell_area`, `max_cell_area`) and the plot settings.
   2. **Test parameters:** check the detection on one image and adjust the parameters if needed.
   3. **Detect & measure:** processes every image.
   4. **Gate:** sets the GFP-positive gate and builds the per-image and per-repeat % positive tables.
   5. **Figures:** swarm plots and histograms (per repeat and whole population), switch response, % positive and fold-change violins, and the superplot.
   6. **Save:** exports every figure (300 dpi png), the tables and the parameters (including the gate) into a time-stamped `results/` folder.

## Acknowledgements

`Functions/func_al_goodplot.m` is adapted from [al_goodplot](https://www.mathworks.com/matlabcentral/fileexchange/91790-al_goodplot-boxblot-violin-plot) on MATLAB Central File Exchange. The statistics annotations have been commented out.
