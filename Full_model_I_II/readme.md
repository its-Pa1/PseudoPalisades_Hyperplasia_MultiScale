# Full Model I & II - Readme

## Main Files/Scripts (for Figures 3.2 and 3.3 of the Thesis)

### Main Scripts
- **main_for_plot_set_6_i.m** (i=F,H): These are the main MATLAB scripts that call all the functions and simulate the 2-D model (3.5.1 of the thesis) on a square grid. They save the results in the folders `Plots_set_6i` and `Videos_set_6i`.
  - **6F** corresponds to experiment 3.1 (Figure 3.2).
  - **6H** corresponds to experiment 3.2 (Figure 3.3).

## Functions
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