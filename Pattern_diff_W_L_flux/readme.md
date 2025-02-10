# Pattern Difference with Limited Flux - Readme

## Main Files/Scripts (for Figure 3.6: 3rd to 6th row of the Thesis)

### Main Scripts
- **main_diff_W_L_flux_i.m** (i=const, I): These are the main MATLAB scripts that call all the functions and simulate the dimensional 1-D model versions. They compute the difference between the full model (3.5.1) and the model version where eq. 3.5.1a is replaced by eq. 3.5.6 (i.e., with gamma_2=0 and the denominator of ph-taxis =1). The results are saved in the folder `Plots_pattern_W_L_flux`. Additionally, they save the patterns for both the full model and the model without self-diffusion.

### Subscripts
- **const**: Corresponds to constant D_T.
- **I**: Corresponds to D_T given by eq. 2.5.15.

## Functions
- **compute_full_model_pattern_const**: Computes the full model 3.5.1 in 1D with constant D_T.
- **compute_full_model_pattern_I**: Computes the full model 3.5.1 in 1D with D_T given by eq. 2.5.15 (as shown in Figure 3.6, first row of the 2nd column).
- **compute_W_L_flux_const**: Computes the model version where 3.5.1a is replaced by eq. 3.5.6 (i.e., gamma_2=0 and the denominator of ph-taxis =1) with constant D_T.
- **compute_W_L_flux_I**: Computes the model where 3.5.1a is replaced by eq. 3.5.6 (gamma_2=0 and the denominator of ph-taxis =1) with D_T given by eq. 2.5.15.
- **set_up_const_diff.m**: Assembles the diffusion matrix (implicit part to be applied for the IMEX method) for constant diffusion.

## Methodology
- **Time Discretization**: For acidity, ECs, and VEGF equations, the IMEX method is used. Diffusion is solved implicitly (implicit Euler), while source and taxis/advection terms are solved explicitly (Euler). For the glioma equation, each term is treated explicitly.
- **Space Discretization**: The standard 5-point stencil (central diff) is used for diffusion terms. The upwind scheme (first order) is used for the glioma and EC taxis terms.

## Requirements
- **Software**: MATLAB (2020a or later version, it should work for older versions as well)