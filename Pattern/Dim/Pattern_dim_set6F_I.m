clear all;
clc;
close all;
% case 6F dimensional
tic
%% parameters in the model(dim form)

% Parameters related to Glioma cells
scale = 60*60*24*30; % scale month
s = (15/3600)*scale; % speed of glioma cells(five times smaller than maximum speed)
lambda = 0.1*scale; % glioma cells turning rate
DT_scale = ((s^2)/lambda); % s^2/lambda
K_M = 0.8; % glioma carrying capacity
mu1 = (0.2/(60*60*24))*scale; % glioma profileration rate due to acidity
mu2 = (0.95/(60*60*24))*scale; % % glioma profileration rate due to VEFG
alpha = 1e+3; %advection constant
gamma_1 = 1e-4; % coefficient of taxis away from acidity
gamma_2 = 1e-4; % coefficient of non-linear diffusion term
%%

% Parameters related to acidity
zeta_S = 1e-7*scale; %acidity updtake rate
gamma_S  = 1e-9*scale; %acidity production rate
D_S =  50*DT_scale; %acidity diffusion i.e 100 times higher than tumor diff
K_S = 10^(-6.4); %maximum acidity

XX = sqrt(D_S/mu1);
C1 = (K_S/XX)^2; % summand with |\nabla h|^2
C2 = (K_M/XX)^2; % summand with |\nabla M|^2

% Parameters related to EC
sigma = (20/3600)*scale ;% speed of endothelial cells 15\mu m per hour
eta = 0.01*scale; % endothelial cells turning rate(10 times faster than glioma ) !!!
D_W = (sigma^2)/(2*eta); % diff coeff of EC
%%%mu_W = 0.56e-6*scale; % proliferation rate of endothelial cells
mu_W = (0.03/(60*60*24))*scale;
chi_W = 0.02; %tactic sensivity of EC towards glioma cells !!!
K_W = 0.3; % EC carrying capacity !!!

% Parameters related to VEGF
K_G = 0.8e-7; %VEGF carrying capacity
D_G  = scale;% 10*scale; % diff coeff of VEGF
gamma_G = 1e-10*scale;%1e-24*scale; % VEGF production rate
zeta_G = 1e-10*scale;%1e-6*scale; % VEGF uptake rate

%% Grid generation
x_end = 1000;
t_end = 15;
dx = 5;
x = 0:dx:1000;
Lx = length(x);
dt = 0.001;
t  = 0:dt:t_end;

%% Computation of DT
%DT = 0.5*DT_scale; % glioma diffusion coeff.
a_d = 1; % aplha for the d(x)
DT = zeros(1,length(x));
for i = 1:length(x)
    if ((250<=x(i)) && (x(i)<=500))
        DT(i) = (sin(0.04*pi*(x(i)-250)))^(2*a_d);
    elseif ((750<=x(i)) && (x(i)<=1000))
        DT(i) = (sin(0.04*pi*(x(i)-750)))^(2*a_d);
    end
end
DT = DT_scale*DT*0.5;
%%  Initial conditions
centerx = 500;

% for tumor
sigma1  = 30;
exponent = ((x-centerx).^2)./(2*sigma1^2);
% M_old       =  K_M*0.1*(exp(-exponent));
M_old       =  0.005*(exp(-exponent));

% for acidity
sigma1  = 20;
exponent = ((x-centerx).^2)./(2*sigma1^2);
S_old       =  (10^(-7))*exp(-exponent);

% initial EC
exponent11 = ((x-950).^2 )./(0.002);
exponent22 = ((x-50).^2 )./(0.002);
% W_old = 0.02*K_W*(exp(-exponent11)+exp(-exponent22));
W_old = 0.0025*(exp(-exponent11)+exp(-exponent22));

% Initial growth factor
exponent = ((x-centerx).^2)./(2*sigma1^2);
G_old       = (10^(-12)*exp(-exponent));

%% Memory allocations
M_total = zeros(Lx,length(t)); % saves M at each time step
S_total = zeros(Lx,length(t)); % saves S at each time step
W_total = zeros(Lx,length(t)); % saves W at each time step
G_total = zeros(Lx,length(t)); % saves G at each time step
M_new = zeros(Lx,1);
flux_M = zeros(Lx,1);
flux_W = zeros(Lx,1);
total_diff_M = zeros(Lx,1);
%% setting up diffusion matrix for const diff to use IMEX method
A_S = set_up_const_diff(x,dt,D_S);
A_W = set_up_const_diff(x,dt,D_W);
A_G = set_up_const_diff(x,dt,D_G);
B_S = zeros(1,size(A_S,2)); % RHS for acidity
B_W = zeros(1,size(A_W,2)); % RHS for ECs
B_G = zeros(1,size(A_G,2)); % RHS for VEGF

%% save a copy of initial condition
M_total(:,1) = M_old;
S_total(:,1) = S_old;
W_total(:,1) = W_old;
G_total(:,1) = G_old;

%% Time loop
for j = 2:length(t)
    
    % glioma diff
    for i = 2:Lx-1
        total_diff_M(i) = DT(i) + ((((DT_scale) - DT(i))*gamma_2*alpha*M_old(i))*...
            (1/(sqrt(C2 + (((M_old(i+1)-M_old(i-1))/(2*dx)).^2)))));
        
        %glioma advective flux
        flux_M(i) = alpha*((DT_scale) - DT(i))*gamma_1*((S_old(i+1)-S_old(i-1))/(2*dx))/...
            ((sqrt(C1 + (((S_old(i+1)-S_old(i-1))/(2*dx))^2))));
        %EC advective flux
        flux_W(i) = -D_W*eta*chi_W* (K_G/ ( (K_G + G_old(i))^2))* ( G_old(i+1)-G_old(i-1))/(2*dx);
    end
    
    
    total_diff_M(1) = DT(1) + (((DT_scale) - DT(1))*gamma_2*alpha*M_old(1));
    total_diff_M(end) = DT(end) + (((DT_scale) - DT(end))*gamma_2*alpha*M_old(end));
    
    flux_M(2) = flux_M(1);
    flux_M(end) = flux_M(end-1);
    flux_W(2) = flux_W(1);
    flux_W(end) = flux_W(end-1);
    
    for i = 2:Lx-1
        % entire glioma eq. using first order upwind for taxis
        if (flux_M(i)<=0)
            M_new(i) = M_old(i) + dt* ( ((0.5/(dx*dx))* ((total_diff_M(i)+ total_diff_M(i+1))*(M_old(i+1)-M_old(i))...
                -(total_diff_M(i-1)+total_diff_M(i))*(M_old(i)-M_old(i-1))))...% diffusion
                + ((1/dx)*(flux_M(i)*M_old(i)-flux_M(i-1)*M_old(i-1)))... % uwpind taxis
                + M_old(i)*(1 - (M_old(i)/K_M))*(mu1*( 1- (S_old(i)/K_S))) ...
                + (((mu2*M_old(i)*W_old(i)/K_W))/( 1 + (W_old(i)/K_W)))); % source
        else
            M_new(i) = M_old(i) + dt* ( ((0.5/(dx*dx))* ((total_diff_M(i)+ total_diff_M(i+1))*(M_old(i+1)-M_old(i))...
                -(total_diff_M(i-1)+total_diff_M(i))*(M_old(i)-M_old(i-1))))... % diffusion
                + ((1/dx)*(flux_M(i+1)*M_old(i+1)-flux_M(i)*M_old(i)))... % uwpind taxis
                + M_old(i)*(1 - (M_old(i)/K_M))*(mu1*( 1- (S_old(i)/K_S))) ...
                + (((mu2*M_old(i)*W_old(i)/K_W))/( 1 + (W_old(i)/K_W)))); % source
        end
        
        if (flux_W(i)<=0)
            B_W(i-1) = (W_old(i)/dt) +(1/dx)*(flux_W(i)*W_old(i)-flux_W(i-1)*W_old(i-1))...% taxis with upwind
                + ((mu_W)*(1-(W_old(i)/K_W))*(G_old(i)/K_G)*W_old(i));%source
        else
            B_W(i-1) = (W_old(i)/dt) + (1/dx)*(flux_W(i+1)*W_old(i+1)-flux_W(i)*W_old(i))... % taxis with upwind
                + ((mu_W)*(1-(W_old(i)/K_W))*(G_old(i)/K_G)*W_old(i));%source
        end
        
        % RHS for acidity & VEGF (source only)
        B_S(i-1) = (S_old(i)/dt) + (gamma_S*M_old(i)/(K_M + M_old(i))) - zeta_S*S_old(i)*(W_old(i)/K_W);
        B_G(i-1) = (G_old(i)/dt) + (gamma_G*M_old(i)*(S_old(i)/K_S)/(K_M + M_old(i))) - zeta_G*(W_old(i)/K_W)*G_old(i);
        
        
    end
    sol_S = (A_S)\B_S';
    sol_W = (A_W)\B_W';
    sol_G = (A_G)\B_G';
    
    % boundary conditions
    M_new(1) = M_new(2);
    M_new(end) = M_new(end-1);
    M_total(:,j) = M_new;
    
    S_total(2:end-1,j) = sol_S;
    S_total(1,j) = S_total(2,j);
    S_total(Lx,j) = S_total(Lx-1,j);
    
    W_total(2:end-1,j) = sol_W;
    W_total(1,j) = W_total(2,j);
    W_total(Lx,j) = W_total(Lx-1,j);
    
    G_total(2:end-1,j) = sol_G;
    G_total(1,j) = G_total(2,j);
    G_total(Lx,j) = G_total(Lx-1,j);
    
    % resetting new solution as old for time marching
    M_old = M_total(:,j);
    S_old = S_total(:,j);
    W_old = W_total(:,j);
    G_old = G_total(:,j);
    
end
%% Plots
% plots of glioma, Acidity, ECs and VEGF with time and space

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_pattern_dim'))
    mkdir('Plots_pattern_dim')
end

figure(1)
surf(x,t,M_total')
view(0,90)
colorbar
caxis([0 0.7e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 7e-5)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Glioma cells' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern_dim/Tumor_pattern_dim_6F_I.png');
%%
figure(2)
surf(x,t,S_total')
view(0,90)
colorbar
caxis([0 5e-7])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_dim/Acidity_pattern_dim_6F_I.png');
%%
figure(3)
surf(x,t,W_total')
view(0,90)
colorbar
caxis([0.0001 0.025])
% caxis([0.000001 0.0005])
%caxis([0.3000000001 0.3000003])
% caxis([0.3000000001 0.30000001])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('ECs' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_dim/EC_pattern_dim_6F_I.png');
%%
figure(4)
surf(x,t,G_total')
view(0,90)
colorbar
caxis([4e-8,7e-8])
% caxis([0.000000002 0.000004])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('VEGF' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_dim/VEGF_pattern_dim_6F_I.png');
%%
figure(6)
plot(x,DT, 'Linewidth', 1.5)
xlabel('X' , 'Fontsize', 15);
ylabel('D_T' , 'Fontsize', 15);
axis tight
title('D_T', 'Fontsize', 15);
saveas(gcf,'Plots_pattern_dim/Diff_coeff_I.png');

%%
toc