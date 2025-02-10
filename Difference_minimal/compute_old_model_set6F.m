function [M_total, S_total] = compute_old_model_set6F(x,y,h,t_end, dt)
%This is the function which solves the old model and save the
% glioma and acidity per 3 days in M_total and S_total
%% Grid generation and memory allocations
[X,Y] = meshgrid(x,y); % square grid
Lx = length(x); % length of x and y
N = (Lx-2)^2; % number of unknowns in the domain
flux_x = zeros(Lx,Lx); % flux in x direction
flux_y = zeros(Lx,Lx); % flux in y direction
advection_x = zeros(Lx-2,Lx-2); % advection term x-component
advection_y = zeros(Lx-2,Lx-2);
source_M = zeros(Lx-2,Lx-2);
diff_temp = zeros(Lx-2,Lx-2);
% these three lines to save 10 values of solutions for 1 time interval(month/week/day)
temp = 1/(10*dt);
M_total = zeros(Lx,Lx,t_end*10+1);
S_total = zeros(Lx,Lx,t_end*10+1);
%% data
delta = 0.2; % delta for small q computation, the parameter to combine two distributions
kappa = 3; % it controls the FA of D_w %
center_x_DW = 450; % position of x for D_w cross
center_y_DW = 450; % position of y for D_w cross

% parameters in the model (all are in \sec time scale then scaled to month)
scale = 60*60*24*30; % scale month
alpha = (1e-11)*scale;
beta  = (1e-9)*scale;
s = (15/3600)*scale;
lambda0 = 0.1*scale;
lambda1 = 0.1*scale;
kp = 0.004*scale;
km = 0.01*scale;
DT_scale = ((s^2)/lambda0);
D_s =  50*DT_scale; % i.e 100 times higher than tumor diff
M_max = 0.8;
S_max = 10^(-6.4);
mu0 = (0.2/(60*60*24))*scale;

%% Initial conditions
% for tumor
sigma1  = 25;   sigma2  = 20;   sigma3  = 10;

centerx = 500; centerx2 = 600; centerx3 =300; %centerx3 = 500;
centery = 500; centery2 = 500; centery3 =400; %centery3 = 600;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
M_old       =  M_max*0.1*(exp(-exponent)+exp(-exponent2)+exp(-exponent3));
% M_old       =  0.005*(exp(-exponent)+exp(-exponent2)+exp(-exponent3));
M_old = M_old';
% for acidity
sigma1  = 15;   sigma2  = 10;   sigma3  = 7.5;

exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
S_old       =  10^(-7)*exp(-exponent)+10^(-7)*exp(-exponent2)+10^(-6.4)*exp(-exponent3);
S_old = S_old';

M_total (:,:,1) = M_old;
S_total (:,:,1) = S_old;
%% calling the diffusion matrix functions
[q,Q_value,a,b,c,diff_stencil,FA,DivDT_x,DivDT_y] = set_tumor_diff_explicit(x,y, delta, kappa, ...
    center_x_DW,center_y_DW,DT_scale);
A_s = set_const_diff(x,D_s,dt);
Q = reshape(Q_value(2:end-1,2:end-1),N,1);
%% Time loop
count = 2;
count2 = 1;
t = 0:dt:t_end;

for j = 1:length(t)
    
    for ii = 2:length(x)-1
        for jj = 2:length(y)-1
            % flux calculation in both directions
            flux_x(ii,jj) =  DivDT_x(ii-1,jj-1)+function_g_new(S_old(ii,jj), S_max,lambda0, lambda1,kp,km)*...
                (a(ii,jj)*((S_old(ii,jj)-S_old(ii-1,jj))/h) + ...
                b(ii,jj)*((S_old(ii,jj)-S_old(ii,jj-1))/h));
            
            flux_y(ii,jj) =  DivDT_y(ii-1,jj-1)+function_g_new(S_old(ii,jj), S_max,lambda0, lambda1,kp,km)*...
                (b(ii,jj)*((S_old(ii,jj)-S_old(ii-1,jj))/h) + ...
                c(ii,jj)*((S_old(ii,jj)-S_old(ii,jj-1))/h));
        end
    end
    flux_x(1,:)   = flux_x(2,:);          flux_y(1,:)   = flux_y(2,:);
    flux_x(end,:) = flux_x(end-1,:);      flux_y(end,:) = flux_y(end-1,:);
    flux_x(:,1)   = flux_x(:,2);          flux_y(:,1)   = flux_y(:,2);
    flux_x(:,end) = flux_x(:,end-1);      flux_y(:,end) = flux_y(:,end-1);
    
    
    for kk = 2:Lx-1
        for ll = 2:Lx-1
            
            % explicit 9 stencils
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
            % glioma source
            source_M(kk-1,ll-1) = (mu0)*((1-(M_old(kk,ll)/M_max)).*(1-(S_old(kk,ll)/S_max))).*M_old(kk,ll);
            
        end
    end
    M_sol = diff_temp+dt*(advection_x+advection_y)+ dt*source_M;
    B_M = reshape(M_old(2:end-1,2:end-1),N,1);
    B_S = reshape(S_old(2:end-1,2:end-1),N,1);
    
    
    RHS_S = B_S+ dt*beta*(B_M./(M_max + B_M)) - dt*alpha*B_S;% RHS containing source & uptake for acidity eq
    
    sol2 = A_s\RHS_S; % solution for acidity at interior nodes
    
    % Boundary conditions
    new = [M_sol(1,1:end);M_sol;M_sol(end,1:end)];
    M_new = [new(1:end,1),new, new(1:end,end)];
    
    S_sol = reshape(sol2,Lx-2,Lx-2);
    new_S = [S_sol(1,1:end);S_sol;S_sol(end,1:end)];
    S_new = [new_S(1:end,1),new_S, new_S(1:end,end)];
    
    
    if (mod(count2,temp)==0)
        M_total (:,:,count) = M_old;
        S_total (:,:,count) = S_old;
        count = count+1;
        
    end
    
    S_old = S_new; % resetting new solution as old for time marching
    M_old = M_new;
    
    count2 = count2+1;
end
