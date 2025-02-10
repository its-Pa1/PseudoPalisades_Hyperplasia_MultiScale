# Difference Self Diffusion - Readme

## Main Files/Scripts (for Figure 5)

### Main Script
- **main_diff_self_d_6F.m**: This is the main MATLAB script that calls all the functions and simulates the 2-D model (in dimensional form) on a square grid. It saves the results (difference between solutions of the full model and the model without flux-limited self-diffusion, gamma_2 = 0, and denominator of ph-taxis = 1) in the folders `Plots_diff_w_sd_6F` and `Videos_diff_w_sd_6F`.

## Functions
- **compute_without_self_diff_6F.m**: Computes the model version 3.5.6 (i.e., with gamma_2 = 0 and denominator of ph-taxis = 1) and saves the involved components: glioma, acidity, EC, and VEGF.
- **compute_new_model_set_6F.m**: Computes the full model 3.5.1 and saves the involved components: glioma, acidity, EC, and VEGF.
- **tissue_Q_macro.m**: Computes the macroscopic tissue Q.
- **function_ag.m**: Computes the a(g) given by eq. 3.3.4.
- **tissue_q_un.m**: Computes the undirected tissue small q given by eq. 2.5.4.
- **set_tumor_data.m**: Stores the tumor-related data like q, DT, FA, div of DT.
- **set_const_diff.m**: Assembles the acidity/EC/Growth factor diffusion matrix to utilize for the IMEX method.
- **set_diff_st.m**: Computes the 9 stencils for anisotropic (and nonlinear) diffusion as mentioned in Table 2.3 of the thesis.

## Methodology
- **Time Discretization**: Explicit method is used for the tumor equation, and the IMEX method is used for other equations. Diffusion is solved implicitly (implicit Euler), while source and taxis/advection terms are solved explicitly (Euler).
- **Space Discretization**: Weikert's 3x3 (central diff, 9 points, given by Table 2.3 of the thesis) stencil discretization (for anisotropic diffusion) is used for tumor diffusion, and the standard 5-point stencil (central diff, given by Table 2.2 of the thesis) is used for acidity/EC/VEGF diffusion.
- **Upwind Scheme**: First-order upwind scheme is used in both x and y directions for the taxis terms (EC and tumor).

## Performance
- **Time**: The code takes approximately 10-15 minutes to run, depending on the machine.

## Requirements
- **Software**: MATLAB (2020a, should work for older versions as well)