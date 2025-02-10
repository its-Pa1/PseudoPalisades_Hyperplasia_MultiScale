Main files/scripts: Figure 3.6: 1st and 2nd row

The file "Pattern_dim_main_6F_i.m" is the main matlab script which calls all the functions and simulates the dimensional 1-D model (eq. 3.5.1) and saves the plots of tumor, acidity, ECs and VEGF w.t.to time and space (for different choices of D_T)

Pattern_dim_main_6F.m: for constant D_T shown in Figure 3.6 (1st column first row)  of the thesis.
Pattern_dim_main_6F_I.m: for D_T(x) (given by eq. 2.5.15) shown in Figure 3.6 (2nd column first row) of the thesis.


Function:
set_up_const_diff.m: assemble the diffusion matrix (implicit part to be applied for IMEX method) for constant diffusion

Method: For acidity, ECs & VEGF equations, IMEX method is used for time discretisation. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler). For glioma equation, each term is treated explicitly. 
 
For space discretisation, standard 5 point stencil(central diff) is used for diffusion terms. Upwind scheme(first order) is used for the glioma and EC taxis terms.


Requirement: Matlab(2020a or later version, it should work for older version also)
