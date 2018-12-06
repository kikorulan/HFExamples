% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex76_3D_40x120x120;

clear all;
close all;

% Dimensions
Nx = 40;
Ny = 120;
Nz = 120;
% Load Initial Pressure
u0_iniMatrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0_ini = matrix2cube(u0_iniMatrix, 240);
% resize
u0 = imresizen(u0_ini, 1/2);
u0Matrix = cube2matrix(u0);
dlmwrite('input_data/initial_pressure_veins_40x120x120.dat', u0Matrix, 'delimiter', ' ');
plot_projection(u0, 1);

% Sound speed
c0 = 1580.0001;
c = c0*ones(Nx, Ny, Nz);
c_matrix = cube2matrix(c);
dlmwrite('input_data/sound_speed.dat', c_matrix, 'delimiter', ' ');
