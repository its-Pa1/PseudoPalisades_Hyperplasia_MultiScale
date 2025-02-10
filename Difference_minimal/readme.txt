Main files/scripts: (for Figure 3.4 of the thesis)

For Figure 3.4: 1st and 2nd column:
The file "main_diff_minimal_set_6F.m" is the main matlab script which calls all the functions and simulates the 2-D model(both new (eq. 3.5.1) and old (eq. 3.5.3) in dimensional form) on a square grid and saves the results(difference between solutions of new model and old model: tumor_new_model-tumor_old_model, acidity_new_model-acidity_old_model) in the folder Plots_diff_minimal_set_6F and Videos_diff_minimal_set_6F. 

For Figure 3.4: 3rd and 4th column:
The file "main_diff_minimal_set_6F_II.m" is the main file to compute the difference between full model (3.5.1) and the system given in the equation 3.5.4 of the thesis and saves the plots and videos in the folder Plots_diff_minimal_set_6F_II and Videos_diff_minimal_set_6F_II


Functions:
compute_old_model_set6F.m: computes the model version 3.5.3 and saves the involved components: glioma and acidity 
compute_full_model_set6F.m: computes the full model 3.5.1 and saves the involved components: glioma, acidity, EC and VEGF
compute_minimal_model_set6F.m : computes the model given in equation 3.5.4 of the thesis and saves the involved components:glioma, acidity.

tissue_Q_macro.m: computes the macroscopic tissue Q
function_ag.m: computes the a(g) given by eq. 3.3.4.
tissue_q_un.m: computes the undirected tissue small q gievn by eq. 2.5.4
set_tumor_data.m: stores the tumor related data like q,DT,FA, div of DT
set_const_diff.m: assemble the acidity/EC/Growth factor diffusion matrix to utilise for IMEX method
set_diff_st.m : computes the 9 stencils for anisotropic(and nonlinear) diffusion as mentioned in Table 2.3 of the thesis.
function_g_new: computes the function g(S) involved in the old model version 3.5.3.
set_tumor_diff_explicit: computes the diff stencils element for explicit Euler method. This function is similar to the "set_diff_st.m" however for a linear case.

Method: for time discretisation explicit method is used for tumor equation and  IMEX method is used for other equations. Diffusion is solved implicitly(implicit Euler) while source and taxis/advection terms explicitly(Euler).
 
For space discretisation, Weikert's 3X3 (central diff, 9 points, given by Table 2.3 of the thesis) stencil discretisation(for anisotropic diffusion) has been used for tumor diffusion and standard 5 point stencil(again central diff, given by Table 2.2 of the thesis)  for acidity/EC/VEGF diffusion.

Upwind scheme(first order) is used in both x and y direction for the taxis terms (EC and tumor).

Time: the codes takes 10-15 minutes depending on the machine


Requirement: Matlab(2020a, it should work for older version also)
