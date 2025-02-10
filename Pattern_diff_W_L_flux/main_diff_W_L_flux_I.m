clear all;
clc;
close all;
%%
dx = 5;
dt = 0.001;
x = 0:dx:1000;
t = 0:dt:15;
%%
[DT,M_full, S_full, W_full, G_full] = compute_full_model_pattern_I();
[M_WL, S_WL, W_WL, G_WL] = compute_W_L_flux_I();
%%
M_diff = M_full - M_WL;
S_diff = S_full - S_WL;
W_diff = W_full - W_WL;
G_diff = G_full - G_WL;

%% Plots
% plots of glioma, Acidity, ECs and VEGF with time and space

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_pattern_W_L_flux'))
    mkdir('Plots_pattern_W_L_flux')
end

figure(1)
surf(x,t,M_WL')
view(0,90)
colorbar
caxis([0 0.7e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 7e-5)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time ' , 'Fontsize', 15);
title('Glioma cells without limited flux' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern_W_L_flux/Tumor_pattern_W_L_flux_I.png');
%%
figure(2)
surf(x,t,S_WL')
view(0,90)
colorbar
caxis([0 5e-7])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity without limited flux' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/Acidity_pattern_W_L_flux_I.png');
%%
figure(3)
surf(x,t,W_WL')
view(0,90)
colorbar
caxis([0.0001 0.025])
%caxis([0.3000000001 0.3000001])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('ECs without limited flux' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/EC_pattern_W_L_flux_I.png');
%%
figure(4)
surf(x,t,G_WL')
view(0,90)
colorbar
caxis([4e-8,7e-8])
% caxis([1.6e-8, 1.5e-6])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('VEGF without limited flux' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/VEGF_pattern_W_L_flux_I.png');
%%
figure(5)
surf(x,t,M_diff')
view(0,90)
colorbar
caxis([-2e-5 2e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 2e-4)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time ' , 'Fontsize', 15);
title('Glioma cells difference' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern_W_L_flux/Tumor_pattern_diff_I.png');
%%
figure(6)
surf(x,t,S_diff')
view(0,90)
colorbar
caxis([-5e-8 5e-8])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity difference' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/Acidity_pattern_diff_I.png');
%%
figure(7)
surf(x,t,W_diff')
view(0,90)
colorbar
caxis([-0.0001 0.01])
%caxis([0.3000000001 0.3000001])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('ECs difference' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/EC_pattern_diff_I.png');
%%
figure(8)
surf(x,t,G_diff')
view(0,90)
colorbar
caxis([-1e-10,1e-8])
% caxis([1.6e-8, 1.5e-6])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('VEGF difference' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/VEGF_pattern_diff_I.png');
%%
figure(9)
plot(x,DT, 'Linewidth', 1.5)
xlabel('X' , 'Fontsize', 15);
ylabel('D_T' , 'Fontsize', 15);
axis tight
title('D_T', 'Fontsize', 15);
saveas(gcf,'Plots_pattern_W_L_flux/Diff_coeff_I.png');

%%
figure(10)
surf(x,t,M_full')
view(0,90)
colorbar
caxis([0 0.7e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 2e-4)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time ' , 'Fontsize', 15);
title('Glioma cells full model' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern_W_L_flux/Tumor_pattern_full_model_I.png');
%%
figure(11)
surf(x,t,S_full')
view(0,90)
colorbar
caxis([0 5e-7])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity full model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/Acidity_pattern_full_model_I.png');
%%
figure(12)
surf(x,t,W_full')
view(0,90)
colorbar
caxis([0.0001 0.025])
%caxis([0.3000000001 0.3000001])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('ECs full model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/EC_pattern_full_model_I.png');
%%
figure(13)
surf(x,t,G_full')
view(0,90)
colorbar
caxis([4e-8,7e-8])
% caxis([1.6e-8, 1.5e-6])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('VEGF full model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_W_L_flux/VEGF_pattern_full_model_I.png');
