clear all;
clc;
close all;
tic
%% space and time related data
h = 5; % spatial spacing
x = 0:h:1000; % x domain
y = x; % y,same as x for a square domain
t_end = 15; % end time %
dt = 0.01;

% old and new model computation
[M_total_full,S_total_full,W_total_full,G_total_full]= compute_full_model_set6F(x,y,h,t_end,dt); % glioma, acidity, EC,VEGF,Q from new model
[M_total_without_sd,S_total_without_sd,W_total_without_sd,G_total_without_sd]= compute_without_self_diff_6F(x,y,h,t_end,dt);
%%
toc % computational time
tic
%% The difference between solutions of new and old model (New-Old)
M_diff3 = M_total_full - M_total_without_sd;
S_diff3 = S_total_full - S_total_without_sd;
W_diff3 = W_total_full - W_total_without_sd;
G_diff3 = G_total_full - G_total_without_sd;
%% videos
% These two loops save the solution video per each 3 days

% In this section Videos are saved in Videos_diff_w_sd_6F folder
% To check if Videos_diff_w_sd_6F folder exists otherwise to create
if not(isfolder('Videos_diff_w_sd_6F'))
    mkdir('Videos_diff_w_sd_6F')
end

% to check which video profile supports available in the machine
% if mp4 is not supported then avi format will be used
profiles = VideoWriter.getProfiles();
check_mp4_support = find(ismember({profiles.Name},'MPEG-4'));

if isempty(check_mp4_support)
    video_ext = '.avi';
    v_pro  = 'Motion JPEG AVI';
else
    video_ext = '.mp4';
    v_pro = 'MPEG-4';
end

videofile = VideoWriter(strcat('Videos_diff_w_sd_6F/tumor_diff_w_sd_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(M_diff3,3)
    figure(1)
    surf(x,y,M_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 14);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_w_sd_6F/acidity_diff_w_sd_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(S_diff3,3)
    figure(2)
    surf(x,y,S_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 14);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_w_sd_6F/EC_diff_w_sd_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(W_diff3,3)
    figure(4)
    surf(x,y,W_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial ECs difference', 'Fontsize', 14);
    else
        title(['ECs difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_w_sd_6F/VEGF_diff_w_sd_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(G_diff3,3)
    figure(5)
    surf(x,y,G_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial VEGF difference', 'Fontsize', 14);
    else
        title(['VEGF difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);
%%
% to save the tumor evolution from old model with the above parameter set
videofile = VideoWriter(strcat('Videos_diff_w_sd_6F/tumor_without_sd_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(M_total_without_sd,3)
    figure(6)
    surf(x,y,M_total_without_sd(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells (without self diff)', 'Fontsize', 12);
    else
        title(['Glioma cells (without self diff) at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 12);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);
%% This loop saves the results for each 30 days in folder Plots_diff_w_sd_6F in png
% or eps format with days in file names

% to check wheather Plots_diff_w_sd_6F folder exists otherwise it makes a folder Plots_diff_w_sd_6F
if not(isfolder('Plots_diff_w_sd_6F'))
    mkdir('Plots_diff_w_sd_6F')
end

for i = 1:10:size(M_diff3,3)
    figure(9)
    surf(x,y,M_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 14);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_diff_w_sd_6F/Tumordiff_w_sd_6F_%ddays',3*(i-1)),'epsc');%eps
    saveas(gcf,sprintf('Plots_diff_w_sd_6F/Tumordiff_w_sd_6F_%ddays.png',3*(i-1)));%png
    
    figure(10)
    surf(x,y,S_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 14);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_diff_w_sd_6F/Aciditydiff_w_sd_6F_%ddays',3*(i-1)),'epsc');
    saveas(gcf,sprintf('Plots_diff_w_sd_6F/Aciditydiff_w_sd_6F_%ddays.png',3*(i-1)))
    
    
    
    figure(11)
    surf(x,y,W_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial ECs difference', 'Fontsize', 14);
    else
        title(['ECs difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_diff_w_sd_6F/Tumordiff_w_sd_6F_%ddays',3*(i-1)),'epsc');%eps
    saveas(gcf,sprintf('Plots_diff_w_sd_6F/ECdiff_w_sd_6F_%ddays.png',3*(i-1)));%png
    
    figure(12)
    surf(x,y,G_diff3(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial VEGF difference', 'Fontsize', 14);
    else
        title(['VEGF difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots_diff_w_sd_6F/Aciditydiff_w_sd_6F_%ddays',3*(i-1)),'epsc');
    saveas(gcf,sprintf('Plots_diff_w_sd_6F/VEGFdiff_w_sd_6F_%ddays.png',3*(i-1)))
end
%% uncomment to save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_diff_self_d_6F.mat');
%%
toc % plots and videos saving time