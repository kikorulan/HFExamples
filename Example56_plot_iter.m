% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex56_RT_3D_4x4x4;
%clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', delimiterIn, headerlinesIn);
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
sound_speed_matrix = importdata('input_data/soundSpeed_10p.dat', delimiterIn, headerlinesIn);
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
saveas(gcf, 'output_data/Example56_C.fig');
saveas(gcf, 'output_data/Example56_C', 'png');
%==============================
% Initial Pressure
%==============================
initial_pressure_matrix = importdata('input_data/initialPressure_3balls.dat', delimiterIn, headerlinesIn);
initial_pressure = matrix2cube(initial_pressure_matrix, Nz);
figure;
surf(x_axis, y_axis, initial_pressure(:, :, 64), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
caxis([0 1]);
box on;
xlabel('y [m]');
ylabel('x [m]');
set(gcf, 'pos', positionY);
saveas(gcf, 'output_data/Example56_U.fig');
saveas(gcf, 'output_data/Example56_U', 'png');

%==================================================
% RECONSTRUCTION 0
%==================================================
recon_matrix = importdata('output_data/pixelPressure_aux.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 10, dx);

recon_matrix = importdata('output_data/zBar_k.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
plot_pixel(reconRT, 15, dx);

%==================================================
% RECONSTRUCTION 1
%==================================================
recon_matrix = importdata('output_data/pixelPressure_56.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 2
%==================================================
recon_matrix = importdata('output_data/pixelPressure_112.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 3
%==================================================
recon_matrix = importdata('output_data/pixelPressure_168.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 4
%==================================================
recon_matrix = importdata('output_data/pixelPressure_504.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)


%==========================================================================================
% VIDEO
%==========================================================================================
% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;
[h, h1, h2] = plot_pixel(initial_pressure, dx, false);
saveas(gcf, 'output_data/Example56_U_3D', 'png');
saveas(gcf, 'output_data/Example56_U_3D.fig');
%==================================================
% 3D render
%==================================================
[h, h1, h2] = plot_pixel(initial_pressure, dx, false, reconRT);
saveas(gcf, 'output_data/Example56_RTrecon_3D', 'png');
saveas(gcf, 'output_data/Example56_RTrecon_3D.fig');
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



 

