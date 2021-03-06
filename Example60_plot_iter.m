% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex60_3D_4balls;
%clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%==========================================================================================
% FORWARD PROBLEM
%==========================================================================================
% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;
positionY = [700 700 600 500];
%==============================
% SOUND SPEED
%==============================
sound_speed_matrix = importdata('input_data/sound_speed.dat', delimiterIn, headerlinesIn);
sound_speed = matrix2cube(sound_speed_matrix, Nz);

figure;
surf(x_axis, y_axis, sound_speed(:, :, 64), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y [m]');
ylabel('x [m]');
box on;
set(gcf, 'pos', positionY);
saveas(gcf, 'output_data/Example60_C.fig');
saveas(gcf, 'output_data/Example60_C', 'png');

%==============================
% Initial Pressure
%==============================
initial_pressure_matrix = importdata('input_data/initial_pressure_4balls.dat', delimiterIn, headerlinesIn);
initial_pressure = matrix2cube(initial_pressure_matrix, Nz);
plot_pixel(initial_pressure, 10, dx);

%==================================================
% RECONSTRUCTION 0
%==================================================
recon_matrix = importdata('output_data/pixelPressure_FISTA_50iter.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
sortCube = sort(reconRT(:));
maxCube = sortCube(end-128);
reconRT = reconRT/maxCube;
plot_pixel(reconRT, 5, dx);
view(41, 6);
ax = gca;
ax.GridAlpha = 0.5;
grid on;
saveas(gcf, 'Example60_RT_recon_FISTA_50iter.fig');
saveas(gcf, 'Example60_RT_recon_FISTA_50iter', 'png');




recon_matrix = importdata('output_data/01.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/02.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/03.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/04.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/05.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/06.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/07.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/08.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/09.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

recon_matrix = importdata('output_data/10.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 5, dx);

%%  para.maxIter = 8000;
%%  para.constraint = 'positivity';
%%  para.outputFL = true;
%%  reconRT2 = conTVdenoising(reconRT, 1e-3, para);



%==================================================
% RECONSTRUCTION AUX
%==================================================
recon_matrix = importdata('output_data/pixelPressure_560.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
sortCube = sort(reconRT(:));
maxCube = sortCube(end-128);
reconRT = reconRT/maxCube;
plot_pixel(reconRT, 5, dx);
view(41, 6);
ax = gca;
ax.GridAlpha = 0.5;
grid on;
saveas(gcf, 'Example60_RT_recon_SPDHG_10epochs.fig');
saveas(gcf, 'Example60_RT_recon_SPDHG_10epochs', 'png');

%==========================================================================================
% VIDEO
%==========================================================================================
% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;

%==================================================
% 3D render
%==================================================
[h, h1, h2] = plot_pixel(initial_pressure, dx, false, reconRT);
saveas(gcf, 'output_data/Example60_RTrecon_3D', 'png');
saveas(gcf, 'output_data/Example60_RTrecon_3D.fig');
% Create New Video with the Image Sequence
outputVideo = VideoWriter('Reconstruction_3D.avi');
outputVideo.FrameRate = 3;
open(outputVideo);

% Loop through the image sequence, load each image, and then write it to the video.
for i = 60:1:120
    disp(i);
    view(h1, i, 35);
    view(h2, i, 35);
    saveas(h, 'FrameVideo', 'jpg');
    img = imread('FrameVideo.jpg');
    writeVideo(outputVideo,img);
end
close(outputVideo);

%==================================================
% 2D render
%==================================================
% Create New Video with the Image Sequence
outputVideo = VideoWriter('Reconstruction_2D.avi');
outputVideo.FrameRate = 10;
open(outputVideo);

% Loop through the image sequence, load each image, and then write it to the video.
ip_norm = initial_pressure(:, :, :)/max(initial_pressure(:));
recon_norm = reconRT(:, :, :)/max(reconRT(:));
% Plot figure
h = figure;
h.Position(3:4) = [1300 500];
h.PaperUnits = 'points';
h.PaperPosition(3:4) = [1300 500];
h1 = subplot(1, 2, 1);
surf(x_axis, y_axis, zeros(Nx, Ny), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
caxis([0 1]);
colorbar();
xlabel('y [m]');
ylabel('x [m]');
h2 = subplot(1, 2, 2);
surf(x_axis, y_axis, zeros(Nx, Ny), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
caxis([0 1]);
colorbar();
xlabel('y [m]');
ylabel('x [m]');


for i = 1:128
    disp(i);
    ip_layer = ip_norm(:, :, i);    
    recon_layer = recon_norm(:, :, i); 
    h1.Children.ZData = ip_layer;
    h2.Children.ZData = recon_layer;
    saveas(h, 'FrameVideo', 'jpg');
    img = imread('FrameVideo.jpg');
    writeVideo(outputVideo,img);
end
close(outputVideo);



 

