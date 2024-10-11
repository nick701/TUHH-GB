# TUHH-GB: Dynamic Mode Decomposition (DMD) Project

## Project Structure

### Code Files
The `code-files/` directory contains the core MATLAB scripts used in the analysis:

- **DMD.m**: Implements Dynamic Mode Decomposition (DMD) for time-series data.
- **DMD_prediction_error.m**: Computes the prediction error using DMD results.
- **DMD_window_size_NRMSE.m**: Analyzes the impact of window size on the Normalized Root Mean Square Error (NRMSE).
- **plots.m**: Generates visualizations of the DMD results.
- **video_plots.m**: Produces video visualizations for dynamic simulations.

### Data Files
The `mat-files/` directory includes `.mat` files that store the datasets and results from the DMD analysis:

- **Seq1.mat**, **Seq2.mat**, **Seq3.mat**: Time-series datasets used for DMD analysis.
- **NRMSE_*.mat**: NRMSE prediction error results for various sequences.
- **f_dmd_standard_*.mat**, **f_dmd_aug_*.mat**: Results from different DMD methods, including standard and augmented DMD.
- **time.mat**: Time data associated with the sequences.