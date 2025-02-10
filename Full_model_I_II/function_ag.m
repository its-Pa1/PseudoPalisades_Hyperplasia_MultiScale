function y = function_ag(chi_w,M,K_M,g,K_G)

% This function computes the coefficient function a(G) given by eq. 3.3.4
% of the thesis.

y = chi_w*(K_G/((K_G+g)^2));

end