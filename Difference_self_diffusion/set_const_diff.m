function A_s = set_const_diff(x,D_s,dt)
% This function computes the diffusion matrix for constant diffusion
% e.g acidity(S)/Endothelial cells/VEGF, here usual 5 point
% stencils has been used for central diff discretization
% D_s =  acidity/Endothelial cells//VEGF diffusion coefficient

h = x(2)-x(1);
Lx = length(x);
N = (Lx-2)^2;

gamma = D_s*dt/(h*h);% the elements of the matrix

% upper most and lower most diagonal
temp = gamma*ones(N,1);

%sub- diagonal
temp1 = gamma*ones(N,1);
temp1(Lx-2:Lx-2:end)=0;

% super-diagonal
temp2 = gamma*ones(N,1);
temp2(Lx-1:Lx-2:end)=0;

% main diagonal
diag_s = (1+4*gamma)*ones(N,1);
diag_s(1:Lx-2) = 1+3*gamma;
diag_s(1) = 1+2*gamma;
diag_s(Lx-2) = 1+2*gamma;
diag_s(1+Lx-2:Lx-2:N-2*(Lx-2)+1)=1+3*gamma;
diag_s(N-(Lx-2)+1) = 1+2*gamma;
diag_s(N-(Lx-2)+2:N-1)=1+3*gamma;
diag_s(2*(Lx-2):Lx-2:(Lx-2-1)*(Lx-2))=1+3*gamma;
diag_s(N) = 1+2*gamma;

% the matrix
A_s = spdiags([-temp, -temp1, diag_s , -temp2, -temp],[-(Lx-2),-1,0,1,Lx-2],N,N);
end