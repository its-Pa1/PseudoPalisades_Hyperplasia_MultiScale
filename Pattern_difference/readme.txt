Main files/scripts: (for Figure 3.7)

The file "Pattern_diff_i_main.m" (i=const, I) is the main Matlab script which calls all the functions and simulates the difference between dimensional 1-D  new model(3.5.1 in the framework of experiment 3.2) and old model (3.5.3) and saves the plots w.t.to time and space in the folder "Plots_pattern_diff". It also saves the patterns for both old and new models.

Subscripts: const, I are the subscripts for different choices of DT(constant and given by eq. 2.5.15 respectively).

Function:
compute_pattern_old_const.m: function computes the old 1D model version 3.5.3 and saves the data for glioma and acidity space-time plots for const DT.

compute_pattern_old_I.m: function computes the old 1D model version 3.5.3 and saves the data for glioma and acidity space-time plots (for degenerate DT given by eq. 2.5.15 and as shown in Figure 3.7 right bottom).

compute_pattern_const.m: function to compute the new model (3.5.1, 1D, experiment 3.2) and to save the glioma, acidity, EC and VEGF space-time plots (for const DT).

compute_pattern_new_I.m: function to compute the new model (3.5.1, 1D, experiment 2) and to save the glioma, acidity, EC and VEGF space-time plots (for degenerate DT given by eq. 2.5.15 and as shown in Figure 3.7 right bottom)

set_up_const_diff.m: assemble the diffusion matrix (implicit part to be applied for IMEX method) for constant diffusion for the new model
set_diff_mat_const.m : assemble the diffusion matrix (implicit part to be applied for IMEX method) for constant diffusion for old model.


Method: For acidity, ECs & VEGF equations, IMEX method is used for time discretisation. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler). For glioma equation, each term is treated explicitly. 
 
For space discretisation, standard 5 point stencil(central diff) is used for diffusion terms. Upwind scheme(first order) is used for the glioma and EC taxis terms.


Requirement: Matlab(2020a or later version, it should work for older version also)
