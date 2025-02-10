Main files/scripts: Figure 3.6: 3rd to 6th row

The file "main_diff_W_L_flux_i.m"(i=const,I) is the main matlab script which calls all the functions and simulates the dimensional 1-D model versions. It computes the difference between the full model (3.5.1) and the model version where eq. 3.5.1a is replaced by  eq. 3.5.6 (i.e with gamma_2=0 and the denominator of ph-taxis =1) and saves the plots in the folder "Plots_pattern_W_L_flux". Also it saves the patterns for both the model versions (full and the without self diff model).

Subscripts: const, I are the subscripts for different choices of DTs (constant and given by eq. 2.5.15 respectively).


Function:
compute_full_model_pattern_const: computes the full model 3.5.1 in 1D with const DT.
compute_full_model_pattern_I : computes the full model 3.5.1 in 1D with DT given by eq. 2.5.15 (as shown in Figure 3.6 (first row of 2nd column)).

compute_W_L_flux_const: computes the model version where 3.5.1a is replaced by  eq. 3.5.6 (i.e gamma_2=0 and the denominator of ph-taxis =1) with const DT
compute_W_L_flux_I: computes the model where 3.5.1a is replaced by  eq. 3.5.6(gamma_2=0 and the denominator of ph-taxis =1) with DT given by eq. 2.5.15.

Function:
set_up_const_diff.m: assemble the diffusion matrix (implicit part to be applied for IMEX method) for constant diffusion

Method: For acidity, ECs & VEGF equations, IMEX method is used for time discretisation. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler). For glioma equation, each term is treated explicitly. 
 
For space discretisation, standard 5 point stencil(central diff) is used for diffusion terms. Upwind scheme(first order) is used for the glioma and EC taxis terms.


Requirement: Matlab(2020a or later version, it should work for older version also)
