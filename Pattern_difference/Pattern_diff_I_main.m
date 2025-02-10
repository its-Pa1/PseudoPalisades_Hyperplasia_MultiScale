clear all;
clc;
close all;
%%
dx = 5;
dt = 0.001;
x = 0:dx:1000;
t = 0:dt:15;
%%
[M_old, S_old] = compute_pattern_old_I();
[DT, M_new, S_new, W_new, G_new] = compute_pattern_new_I();
%%
M_diff = M_new - M_old;
S_diff = S_new - S_old;

%% Plots
% plots of glioma, Acidity, ECs and VEGF with time and space

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_pattern_diff'))
    mkdir('Plots_pattern_diff')
end

figure(1)
surf(x,t,M_old')
view(0,90)
colorbar
caxis([0 0.7e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 0.7e-5)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Glioma cells old model' , 'Fontsize', 15);
axis tight
% axis([300 700 0, t(end)]);
saveas(gcf,'Plots_pattern_diff/Tumor_pattern_old_I.png');
%%
figure(2)
surf(x,t,S_old')
view(0,90)
colorbar
caxis([0 5e-7])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity old model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_diff/Acidity_pattern_old_I.png');
%%
figure(3)
surf(x,t,M_new')
view(0,90)
colorbar
caxis([0 0.7e-5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 2e-4)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Glioma cells new model' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern_diff/Tumor_pattern_new_I.png');
%%
figure(4)
surf(x,t,S_new')
view(0,90)
colorbar
caxis([0 5e-7])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity new model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_diff/Acidity_pattern_new_I.png');
%%
figure(5)
surf(x,t,W_new')
view(0,90)
colorbar
caxis([0.0001 0.025])
% caxis([0.300000000001 0.300000001])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('ECs new model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_diff/EC_pattern_new_I.png');
%%
figure(6)
surf(x,t,G_new')
view(0,90)
colorbar
caxis([4e-8,7e-8])
% caxis([1.6e-8, 1.5e-6])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('VEGF new model' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_diff/VEGF_pattern_new_I.png');
%%
figure(7)
surf(x,t,M_diff')
view(0,90)
colorbar
caxis([-2e-6 2e-6]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 2e-4)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Glioma cells difference' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern_diff/Tumor_pattern_diff_I.png');
%%
figure(8)
surf(x,t,S_diff')
view(0,90)
colorbar
caxis([-1e-7 1e-7])
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity difference' , 'Fontsize', 15);
axis tight
saveas(gcf,'Plots_pattern_diff/Acidity_pattern_diff_I.png');
%%
figure(9)
plot(x,DT, 'Linewidth', 1.5)
xlabel('X' , 'Fontsize', 15);
ylabel('D_T' , 'Fontsize', 15);
axis tight
title('D_T', 'Fontsize', 15);
saveas(gcf,'Plots_pattern_diff/Diff_coeff_I.png');