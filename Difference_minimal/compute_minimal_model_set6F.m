function [M_total,S_total]= compute_minimal_model_set6F(x,y,h,t_end,dt)
%% Grid generation and memory allocations
[X,Y] = meshgrid(x,y); % square grid
Lx = length(x); % length of x and y
N = (Lx-2)^2; % number of unknowns in the domain
flux_x = zeros(Lx,Lx); % Glioma flux in x direction
flux_y = zeros(Lx,Lx); % Glioma flux in y direction
advection_x = zeros(Lx-2,Lx-2); % Glioma advection term x-component
advection_y = zeros(Lx-2,Lx-2); % Glioma advection term y-component
source_M = zeros(Lx-2,Lx-2); % source glioma
diff_temp = zeros(Lx-2,Lx-2); % temprory diffusion stencils of glioma equation
% these  lines to save 10 values of solutions for 1 time interval(month/week/day)
temp = 1/(10*dt);
M_total = zeros(Lx,Lx,t_end*10+1);
S_total = zeros(Lx,Lx,t_end*10+1);

%% data
delta = 0.2; % delta for small q computation, the parameter to combine two distributions
kappa = 3; % it controls the FA of D_w %
center_x_DW = 450; % position of x for D_w cross
center_y_DW = 450; % position of y for D_w cross

% parameters in the model (all are in \sec time scale then scaled to month)

% Parameters related to Glioma cells
scale = 60*60*24*30; % scale month
s = (15/3600)*scale; % speed of glioma cells
lambda = 0.1*scale; % glioma cells turning rate
DT_scale = ((s^2)/lambda); % s^2/lambda
K_M = 0.8; % glioma carrying capacity %0.7 for 6C
mu1 = (0.2/(60*60*24))*scale; % glioma profileration rate due to acidity 0.2 for 6D
alpha = 1e+3; %advection constant
%%%gamma_1 = 3e-4; % coefficient of taxis away from acidity
gamma_1 = 1e-4; % in 6a war 5e-4; in 6b war 1e-3, M(1-M)
gamma_2 = 1e-4; % coefficient of non-linear diffusion term

% Parameters related to acidity
zeta_S = 1e-11*scale; %acidity uptake rate
gamma_S  = 1e-9*scale; %acidity production rate
D_S =  50*DT_scale; %acidity diffusion i.e 50 times higher than tumor diff
K_S = 10^(-6.4); %maximum acidity

XX = sqrt(D_S/mu1);
C1 = (K_S/XX)^2; % summand with |\nabla h|^2
C2 = (K_M/XX)^2; % summand with |\nabla M|^2

%% Initial conditions
centerx = 500; centerx2 = 600; centerx3 =300; %centerx3 = 500;
centery = 500; centery2 = 500; centery3 =400; %centery3 = 600;

% for tumor
sigma1  = 25;   sigma2  = 20;   sigma3  = 10;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
% M_old       =  0.005*(exp(-exponent)+exp(-exponent2)+exp(-exponent3));
M_old       =  K_M*0.1*(exp(-exponent)+exp(-exponent2)+exp(-exponent3));
M_old = M_old';

% for acidity
sigma1  = 15;   sigma2  = 10;   sigma3  = 7.5;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
S_old       =  10^(-7)*exp(-exponent)+10^(-7)*exp(-exponent2)+10^(-6.4)*exp(-exponent3);
S_old = S_old';


% storing inital conditions
M_total (:,:,1) = M_old;
S_total (:,:,1) = S_old;

%% calling the diffusion matrix functions
[q,Q_old,a,b,c,FA,DivDT_x,DivDT_y] = set_tumor_data(x,y, delta, kappa, ...
    center_x_DW,center_y_DW,DT_scale);
A_S = set_const_diff(x,D_S,dt); % diff matrix for Acidity
%% Time loop
count = 2;
count2 = 1;
t = 0:dt:t_end;

for j = 1:length(t)
    
    for ii = 2:length(x)-1
        for jj = 2:length(y)-1
            
            % glioma advective flux calculation in both directions
            flux_x(ii,jj) =  DivDT_x(ii-1,jj-1)+(alpha*gamma_1/(sqrt(C1+ (((S_old(ii+1,jj)-S_old(ii-1,jj))/(2*h))^2) ...
                + (((S_old(ii,jj+1)-S_old(ii,jj-1))/(2*h))^2))))*...
                ((DT_scale-a(ii,jj))*((S_old(ii+1,jj)-S_old(ii-1,jj))/(2*h)) + ...
                (-b(ii,jj))*((S_old(ii,jj+1)-S_old(ii,jj-1))/(2*h)));
            
            flux_y(ii,jj) =  DivDT_y(ii-1,jj-1)+(alpha*gamma_1/(sqrt(C1+ (((S_old(ii+1,jj)-S_old(ii-1,jj))/(2*h))^2) ...
                + (((S_old(ii,jj+1)-S_old(ii,jj-1))/(2*h))^2))))*...
                (-b(ii,jj)*((S_old(ii+1,jj)-S_old(ii-1,jj))/(2*h)) + ...
                (DT_scale-c(ii,jj))*((S_old(ii,jj+1)-S_old(ii,jj-1))/(2*h)));
            
        end
    end
    % BC for fluxes
    flux_x(1,:)   = flux_x(2,:);          flux_y(1,:)   = flux_y(2,:);
    flux_x(end,:) = flux_x(end-1,:);      flux_y(end,:) = flux_y(end-1,:);
    flux_x(:,1)   = flux_x(:,2);          flux_y(:,1)   = flux_y(:,2);
    flux_x(:,end) = flux_x(:,end-1);      flux_y(:,end) = flux_y(:,end-1);
    
    
    % diffusion stencils
    diff_stencil = set_diff_st(x,a,b,c,M_old,alpha, gamma_2,DT_scale,C2);
    
    for kk = 2:Lx-1
        for ll = 2:Lx-1
            
            % 9 stencils (explicit)
            diff_temp(kk-1,ll-1) = dt*diff_stencil.alpha1(kk-1,ll-1)*M_old(kk-1,ll+1)...
                + dt*diff_stencil.alpha2(kk-1,ll-1)*M_old(kk,ll+1)...
                + dt*diff_stencil.alpha3(kk-1,ll-1)*M_old(kk+1,ll+1)...
                + dt*diff_stencil.alpha4(kk-1,ll-1)*M_old(kk-1,ll)...
                + (1+dt*diff_stencil.alpha5(kk-1,ll-1))*M_old(kk,ll)...
                + dt*diff_stencil.alpha6(kk-1,ll-1)*M_old(kk+1,ll)...
                + dt*diff_stencil.alpha7(kk-1,ll-1)*M_old(kk-1,ll-1)...
                + dt*diff_stencil.alpha8(kk-1,ll-1)*M_old(kk,ll-1)...
                + dt*diff_stencil.alpha9(kk-1,ll-1)*M_old(kk+1,ll-1);
            % upwinding in both directions
            if(flux_x(kk,ll)<=0)
                advection_x(kk-1,ll-1) = (flux_x(kk,ll)*M_old(kk,ll) - flux_x(kk-1,ll)*M_old(kk-1,ll))/h;
            else
                advection_x(kk-1,ll-1) = (flux_x(kk+1,ll)*M_old(kk+1,ll) - flux_x(kk,ll)*M_old(kk,ll))/h;
            end
            
            if (flux_y(kk,ll)<=0)
                advection_y(kk-1,ll-1) = (flux_y(kk,ll)*M_old(kk,ll) - flux_y(kk,ll-1)*M_old(kk,ll-1))/h;
            else
                advection_y(kk-1,ll-1) = (flux_y(kk,ll+1)*M_old(kk,ll+1) - flux_y(kk,ll)*M_old(kk,ll))/h;
            end
            
            source_M(kk-1,ll-1) = M_old(kk,ll)*(1 - (M_old(kk,ll)/K_M))*(mu1*( 1- (S_old(kk,ll)/K_S)));
            
        end
    end
    M_sol = diff_temp + dt*(advection_x + advection_y) + dt*source_M; % glioma cell update (explicit)
    
    B_M = reshape(M_old(2:end-1,2:end-1),N,1);
    B_S = reshape(S_old(2:end-1,2:end-1),N,1);
    
    
    
    RHS_S = B_S + dt*gamma_S*(B_M./(K_M + B_M)) - dt*zeta_S*B_S;% RHS containing source & uptake for acidity eq
    
    sol_S = A_S\RHS_S; % solution for acidity at interior nodes
    
    
    % Boundary conditions
    new = [M_sol(1,1:end);M_sol;M_sol(end,1:end)];
    M_new = [new(1:end,1),new, new(1:end,end)];
    
    S_sol = reshape(sol_S,Lx-2,Lx-2);
    new_S = [S_sol(1,1:end);S_sol;S_sol(end,1:end)];
    S_new = [new_S(1:end,1),new_S, new_S(1:end,end)];
    
    
    if (mod(count2,temp)==0)
        M_total (:,:,count) = M_old;
        S_total (:,:,count) = S_old;
        
        
        count = count+1;
        
    end
    % resetting new solution as old for time marching
    S_old = S_new;
    M_old = M_new;
    
    
    count2 = count2+1;
end

end

