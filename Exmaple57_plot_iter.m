% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex57_RT_3D_5x5x5;
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
% ITERATIVE PROBLEM
%==========================================================================================
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
xlabel('y (m)');
ylabel('x (m)');
title('Sound Speed');

%==============================
% INITIAL PRESSURE
%==============================
initial_pressure_matrix = importdata('input_data/initialPressure_3balls.dat', delimiterIn, headerlinesIn);
initial_pressure = matrix2cube(initial_pressure_matrix, Nz);

%==================================================
% INITIAL RECONSTRUCTION
%==================================================
recon_matrix = importdata('output_data/x0.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
reconRT_pos = reconRT;
reconRT_pos(reconRT < 0) = 0;
plot_cube(reconRT_pos);

%==================================================
% RECONSTRUCTION 0
%==================================================
recon_matrix = importdata('output_data/pixelPressure_aux.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 1
%==================================================
recon_matrix = importdata('output_data/pixelPressure_98.dat', ' ', 0);
reconRT1 = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 2
%==================================================
recon_matrix = importdata('output_data/pixelPressure_196.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 3
%==================================================
recon_matrix = importdata('output_data/pixelPressure_294.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)

%==================================================
% RECONSTRUCTION 4
%==================================================
recon_matrix = importdata('output_data/pixelPressure_392.dat', ' ', 0);
reconRT = matrix2cube(recon_matrix, Nz);
%plot_pixel(reconRT1, dx, false);
plot_cube(reconRT)
