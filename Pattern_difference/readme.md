# Pattern Difference - Readme

## Main Files/Scripts (for Figure 3.7 of the Thesis)

### Main Script
- **Pattern_diff_i_main.m** (i=const, I): This is the main MATLAB script that calls all the functions and simulates the difference between the dimensional 1-D new model (3.5.1 in the framework of experiment 3.2) and the old model (3.5.3). It saves the plots with respect to time and space in the folder `Plots_pattern_diff`. It also saves the patterns for both old and new models.

### Subscripts
- **const**: Corresponds to constant D_T.
- **I**: Corresponds to D_T given by eq. 2.5.15.

## Functions
- **compute_pattern_old_const.m**: Computes the old 1D model version 3.5.3 and saves the data for glioma and acidity space-time plots for constant D_T.
- **compute_pattern_old_I.m**: Computes the old 1D model version 3.5.3 and saves the data for glioma and acidity space-time plots for degenerate D_T given by eq. 2.5.15 (as shown in Figure 3.7, right bottom).
- **compute_pattern_const.m**: Computes the new model (3.5.1, 1D, experiment 3.2) and saves the glioma, acidity, EC, and VEGF space-time plots for constant D_T.
- **compute_pattern_new_I.m**: Computes the new model (3.5.1, 1D, experiment 3.2) and saves the glioma, acidity, EC, and VEGF space-time plots for degenerate D_T given by eq. 2.5.15 (as shown in Figure 3.7, right bottom).
- **set_up_const_diff.m**: Assembles the diffusion matrix (implicit part to be applied for the IMEX method) for constant diffusion for the new model.
- **set_diff_mat_const.m**: Assembles the diffusion matrix (implicit part to be applied for the IMEX method) for constant diffusion for the old model.

## Methodology
- **Time Discretization**: For acidity, ECs, and VEGF equations, the IMEX method is used. Diffusion is solved implicitly (implicit Euler), while source and taxis/advection terms are solved explicitly (Euler). For the glioma equation, each term is treated explicitly.
- **Space Discretization**: The standard 5-point stencil (central diff) is used for diffusion terms. The upwind scheme (first order) is used for the glioma and EC taxis terms.

## Requirements
- **Software**: MATLAB (2020a or later version, it should work for older versions as well)