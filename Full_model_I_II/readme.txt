Main files/scripts:
The file "main_for_plot_set_6_i.m" (i=F,H) are the main matlab script which call all the functions and simulates the 2-D model (3.5.1 of the thesis) on a square grid and saves the results in the folder Plots_set_6i and Videos_set_6i.

subscript 6F corresponds to experiment 3.1, while 6H to experiment 3.2 (for figure 3.2 and 3.3 respectively)

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
