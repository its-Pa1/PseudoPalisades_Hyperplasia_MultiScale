function [M, S] = compute_pattern_old_I()
%% data (dimensionalised parameters)
scale = 60*60*24*30; % scale month
alpha = 1e-11*scale;
beta = 1e-9*scale;
s = (15/3600)*scale;
lambda0 = (0.1)*scale;
lambda1 = 0.1*scale;
kp = 0.004*scale;
km = 0.01*scale;
M_max = 0.8;
S_max = 10^(-6.4);
mu0 = (0.2/(60*60*24))*scale;
DT_scale = ((s^2)/lambda0);
Ds = 50*DT_scale;
%% mesh generation, discretization and memory allocation
h = 5;
x = 0:h:1000;
Lx = length(x);
dt = 0.001;
t = 0:dt:15;
M = zeros(Lx,length(t)); % saves M at each time step
S = zeros(Lx,length(t)); % saves S at each time step
flux = zeros(Lx,1);

%% Setting up both the diff matrices
% tumor diff matrix with 1D DT

[D,A] = set_diff_mat_I(x,dt,DT_scale);

% acidity diff matrix
A2  = diag((1/dt) + (2*Ds/(h*h)) * ones(Lx-2,1), 0) + diag(-Ds/(h*h)* ones(Lx-3,1), -1) ...
    + diag(- (Ds/(h*h))* ones(Lx-3,1), 1);
A2(1,1) = A2(1,1)-Ds/(h*h); % BC
A2(end,end) = A2(end,end)-Ds/(h*h); %BC

B_M = zeros(1,size(A,2)); % RHS for tumor
B_S = zeros(1,size(A,2)); % RHS for acidity
%% Initial conditions
centerx = 500;

% for tumor
sigma1  = 30;
exponent = ((x-centerx).^2)./(2*sigma1^2);
M_old       =  0.005*(exp(-exponent));

% for acidity
sigma1  = 20;
exponent = ((x-centerx).^2)./(2*sigma1^2);
S_old       =  (10^(-7))*exp(-exponent);

%% save a copy of initial condition
M(:,1) = M_old;
S(:,1) = S_old;

%% Time loop
for j= 2:length(t)
    
    % flux
    for i = 2: Lx-1
        flux(i) = (D(i)*function_g_new(S_old(i), S_max,lambda0, lambda1,kp,km)* (S_old(i+1)-S_old(i))*(1/h))...
            +((D(i)-D(i-1))/h);
    end
    flux(1) = flux(2);
    flux(Lx) = flux(Lx-1);
    
    for i = 2:Lx-1
        % computation of RHS for tumor eq(source+ advection)
        if (flux(i)<=0)
            B_M(i-1) = (M_old(i)/dt) +(1/h)*(flux(i)*M_old(i)-flux(i-1)*M_old(i-1))...% taxis with upwind
                + ((mu0)*(1-(M_old(i)/M_max))*(1-(S_old(i)/S_max))*M_old(i));%source
        else
            B_M(i-1) = (M_old(i)/dt) + (1/h)*(flux(i+1)*M_old(i+1)-flux(i)*M_old(i))...
                +( (mu0)*(1-(M_old(i)/M_max))*(1-(S_old(i)/S_max))*M_old(i));
        end
        
        % RHS for acidity (source only)
        B_S(i-1) = (S_old(i)/dt)+(beta*M_old(i)-alpha*S_old(i));
        
    end
    
    % solutions at interior nodes
    sol_M = (A)\B_M';
    sol_S = (A2)\B_S';
    
    % boundary conditions
    M(2:end-1,j) = sol_M;
    S(2:end-1,j) = sol_S;
    M(1,j) = M(2,j);
    M(Lx,j) = M(Lx-1,j);
    S(1,j) = S(2,j);
    S(Lx,j) = S(Lx-1,j);
    
    % resetting new solution as old for time marching
    M_old = M(:,j);
    S_old = S(:,j);
    
end

end