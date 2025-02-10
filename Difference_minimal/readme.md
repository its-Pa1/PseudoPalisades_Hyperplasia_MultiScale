# Difference Minimal - Readme

## Main Files/Scripts (for Figure 3.4 of the Thesis)

### For Figure 3.4: 1st and 2nd Column
- **main_diff_minimal_set_6F.m**: This is the main MATLAB script that calls all the functions and simulates the 2-D model (both new (eq. 3.5.1) and old (eq. 3.5.3) in dimensional form) on a square grid. It saves the results (difference between solutions of new model and old model: tumor_new_model - tumor_old_model, acidity_new_model - acidity_old_model) in the folders `Plots_diff_minimal_set_6F` and `Videos_diff_minimal_set_6F`.

### For Figure 3.4: 3rd and 4th Column
- **main_diff_minimal_set_6F_II.m**: This script computes the difference between the full model (3.5.1) and the system given in equation 3.5.4 of the thesis. It saves the plots and videos in the folders `Plots_diff_minimal_set_6F_II` and `Videos_diff_minimal_set_6F_II`.

## Functions
- **compute_old_model_set6F.m**: Computes the model version 3.5.3 and saves the involved components: glioma and acidity.
- **compute_full_model_set6F.m**: Computes the full model 3.5.1 and saves the involved components: glioma, acidity, EC, and VEGF.
- **compute_minimal_model_set6F.m**: Computes the model given in equation 3.5.4 of the thesis and saves the involved components: glioma and acidity.
- **tissue_Q_macro.m**: Computes the macroscopic tissue Q.
- **function_ag.m**: Computes the a(g) given by eq. 3.3.4.
- **tissue_q_un.m**: Computes the undirected tissue small q given by eq. 2.5.4.
- **set_tumor_data.m**: Stores the tumor-related data like q, DT, FA, div of DT.
- **set_const_diff.m**: Assembles the acidity/EC/Growth factor diffusion matrix to utilize for the IMEX method.
- **set_diff_st.m**: Computes the 9 stencils for anisotropic (and nonlinear) diffusion as mentioned in Table 2.3 of the thesis.
- **function_g_new.m**: Computes the function g(S) involved in the old model version 3.5.3.
- **set_tumor_diff_explicit.m**: Computes the diff stencils element for the explicit Euler method. This function is similar to "set_diff_st.m" but for a linear case.

## Methodology
- **Time Discretization**: Explicit method is used for the tumor equation, and the IMEX method is used for other equations. Diffusion is solved implicitly (implicit Euler), while source and taxis/advection terms are solved explicitly (Euler).
- **Space Discretization**: Weikert's 3x3 (central diff, 9 points, given by Table 2.3 of the thesis) stencil discretization (for anisotropic diffusion) is used for tumor diffusion, and the standard 5-point stencil (central diff, given by Table 2.2 of the thesis) is used for acidity/EC/VEGF diffusion.
- **Upwind Scheme**: First-order upwind scheme is used in both x and y directions for the taxis terms (EC and tumor).

## Performance
- **Time**: The code takes approximately 10-15 minutes to run, depending on the machine.

## Requirements
- **Software**: MATLAB (2020a, should work for older versions as well)