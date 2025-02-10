Main files/scripts: for figure 5

The file "main_diff_self_d_6F.m" is the main matlab script which calls all the functions and simulates the 2-D model(in dimensional form) on a square grid and saves the results(difference between solutions of full model and model without flux limited self diffusion, gamma_2 = 0 and denominator of ph-taxis=1 ) in the folder Plots_diff_w_sd_6F and Videos_diff_w_sd_6F


Functions:
compute_without_self_diff_6F.m: computes the model version 3.5.6(i.e. with gamma_2 = 0 and denominator of ph-taxis = 1)and saves the involved components: glioma, acidity, EC and VEGF. 
compute_new_model_set_6F.m: computes the full model 3.5.1 and saves the involved components glioma, acidity, EC & VEGF.
tissue_Q_macro.m: computes the macroscopic tissue Q
function_ag.m: computes the a(g) given by eq. 3.3.4.
tissue_q_un.m: computes the undirected tissue small q gievn by eq. 2.5.4
set_tumor_data.m: stores the tumor related data like q,DT,FA, div of DT
set_const_diff.m: assemble the acidity/EC/Growth factor diffusion matrix to utilise for IMEX method
set_diff_st.m : computes the 9 stencils for anisotropic(and nonlinear) diffusion as mentioned in Table 2.3 of the thesis.


Method: for time discretisation explicit method is used for tumor equation and  IMEX method is used for other equations. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler).
 
For space discretisation, Weikert's 3X3 (central diff, 9 points, given by Table 2.3 of the thesis) stencil discretisation(for anisotropic diffusion) has been used for tumor diffusion and standard 5 point stencil(again central diff, given by Table 2.2 of the thesis)  for acidity/EC/VEGF diffusion.

Upwind scheme(first order) is used in both x and y direction for the taxis terms (EC and tumor).

Time: the codes takes 10-15 minutes depending on the machine


Requirement: Matlab(2020a, it should work for older version also)