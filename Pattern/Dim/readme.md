# Pattern - Dimensional - Readme

## Main Files/Scripts (for Figure 3.6: 1st and 2nd row of the Thesis)

### Main Scripts
- **Pattern_dim_main_6F.m**: This is the main MATLAB script for constant D_T, shown in Figure 3.6 (1st column, first row) of the thesis. It calls all the functions and simulates the dimensional 1-D model (eq. 3.5.1) and saves the plots of tumor, acidity, ECs, and VEGF with respect to time and space.
- **Pattern_dim_main_6F_I.m**: This script is for D_T(x) (given by eq. 2.5.15), shown in Figure 3.6 (2nd column, first row) of the thesis.

## Functions
- **set_up_const_diff.m**: Assembles the diffusion matrix (implicit part to be applied for the IMEX method) for constant diffusion.

## Methodology
- **Time Discretization**: For acidity, ECs, and VEGF equations, the IMEX method is used. Diffusion is solved implicitly (implicit Euler), while source and taxis/advection terms are solved explicitly (Euler). For the glioma equation, each term is treated explicitly.
- **Space Discretization**: The standard 5-point stencil (central diff) is used for diffusion terms. The upwind scheme (first order) is used for the glioma and EC taxis terms.

## Requirements
- **Software**: MATLAB (2020a or later version, it should work for older versions as well)