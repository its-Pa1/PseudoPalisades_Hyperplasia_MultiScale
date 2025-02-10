clear all;
clc;
close all;
tic
%% space and time related data
h = 5; % spatial spacing
x = 0:h:1000; % x domain
y = x; % y,same as x for a square domain
t_end = 15; % end time %
dt = 0.005;

% old and new model computation
[M_total_old,S_total_old] = compute_old_model_set6F(x,y,h,t_end,dt); % glioma and acidity from old model
[M_total_new,S_total_new,W_total_new,G_total_new]= compute_full_model_set6F(x,y,h,t_end,dt); % glioma, acidity, EC,VEGF,Q from new model
%%
toc % computational time
tic
%% The difference between solutions of new and old model (New-Old)
M_diff1 = M_total_new - M_total_old;
S_diff1 = S_total_new - S_total_old;
%% videos
% These two loops save the solution video per each 3 days

% In this section Videos are saved in Videos_diff_minimal_6F folder
% To check if Videos_diff_minimal_6F folder exists otherwise to create
if not(isfolder('Videos_diff_minimal_6F'))
    mkdir('Videos_diff_minimal_6F')
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

videofile = VideoWriter(strcat('Videos_diff_minimal_6F/tumor_diff_minimal_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(M_diff1,3)
    figure(4)
    surf(x,y,M_diff1(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_minimal_6F/acidity_diff_minimal_6F', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(S_diff1,3)
    figure(5)
    surf(x,y,S_diff1(:,:,i)')
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
%%
% to save the tumor evolution from old model with the above parameter set
videofile = VideoWriter(strcat('Videos_diff_minimal_6F/tumor_old_model', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(M_total_old,3)
    figure(4)
    surf(x,y,M_total_old(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    drawnow
    if (i==1)
        title('Initial glioma cells (minimal)', 'Fontsize', 14);
    else
        title(['Glioma cells (minimal model) at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 14);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);
%% This loop saves the results for each 30 days in folder Plots_diff_minimal_6F in png
% or eps format with days in file names

% to check wheather Plots folder exists otherwise it makes a folder Plots
if not(isfolder('Plots_diff_minimal_6F'))
    mkdir('Plots_diff_minimal_6F')
end


for i = 1:10:size(M_diff1,3)
    figure(9)
    surf(x,y,M_diff1(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots_Videos_diff_minimal_6F/Tumor_diff_minimal_6F_%ddays',3*(i-1)),'epsc');%eps
    saveas(gcf,sprintf('Plots_diff_minimal_6F/Tumor_diff_minimal_6F_%ddays.png',3*(i-1)));%png
    
    figure(10)
    surf(x,y,S_diff1(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots_Videos_diff_minimal_6F/Acidity_diff_minimal_6F_%ddays',3*(i-1)),'epsc');
    saveas(gcf,sprintf('Plots_diff_minimal_6F/Acidity_diff_minimal_6F_%ddays.png',3*(i-1)))
    
end
%% uncomment to save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_diff_minimal_6F.mat');

%%
toc % plots and videos saving time