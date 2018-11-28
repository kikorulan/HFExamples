% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex72_3D_veins_heterogeneous;

clear all;
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

%============================================================================================================================================
% FORWARD PROBLEM
%============================================================================================================================================


% SOUND SPEED

% Import data
sound_speed = importdata('output_data/sound_speed.dat', ' ', 0);
sound_speed_cube = matrix2cube(sound_speed, Nz);
plot_cube(sound_speed_cube);

%==================================================
% TIME SIGNAL - RT
%==================================================
% Import data
filenameData = 'output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

% Import data
filenameData = 'input_data/forwardSignal_reference_1600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

%==================================================
% TIME SIGNAL - kWave
%==================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex68_3D_veins_resize;
sensor_data_hom = h5read('output_data/Example68_forward_output_1600sensors.h5', '/p');
inputKWave = sensor_data_hom(:, :);
% Plot
figure;
imagesc(inputKWave);
box on;
colorbar();

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex72_3D_veins_heterogeneous;
sensor_data_het = h5read('output_data/Example72_forward_output_1600sensors.h5', '/p');
inputKWave = sensor_data_het(:, :);
% Plot
figure;
imagesc(inputKWave);
box on;
colorbar();

%==================================================
% Compute forward signal difference
%==================================================
% Import data
filenameData = 'input_data/forwardSignal_1600sensors_dt1.6e-8_lowres_interp_delay2.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
signalRT_1 = timeSignal(2:end, :);
% Import data
filenameData = 'input_data/forwardSignal_kWave_adjoint_1600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
signalRT_2 = timeSignal(2:end, :);
% Create signal
differenceRT = signalRT_1 - signalRT_2;
forwardDif = [timeRT; differenceRT];
dlmwrite('input_data/forwardSignal_difference_1600sensors.dat', forwardDif, 'delimiter', ' ');

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================
% Position
position     = [700 700 300 630];
positionY    = [700 700 320 630];
positionBar  = [700 700 363 630];
positionYBar = [700 700 390 630];
set(0,'DefaultFigurePaperPositionMode','auto');

%==================================================
% INITIAL PRESSURE
%==================================================
% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
saveas(gcf, 'figures/Example68_initial_pressure.fig');
% Load Initial Pressure
load ./input_data/initial_pressure_veins_smooth;
%plot_projection(initial_pressure_veins_smooth, dx);

%==================================================
% Reconstruction - kWave
%==================================================
% 57600 sensors
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex68_3D_veins_resize;
p0_recon_PML = h5read('output_data/Example68_adjoint_output_57600sensors.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
p0_recon_hom = p0_recon/max(p0_recon(:));

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex72_3D_veins_heterogeneous;
p0_recon_PML = h5read('output_data/Example72_adjoint_output_57600sensors.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
p0_recon_het = p0_recon/max(p0_recon(:));
plot_projection(p0_recon_het-p0_recon_hom, 1);

% 1600 sensors
p0_recon_PML = h5read('output_data/Example72_adjoint_output_1600sensors.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
p0_recon = p0_recon/max(p0_recon(:));
plot_projection(p0_recon, dx);

%==================================================
% Reconstruction - RT
%==================================================
% Import data
pixelPressureMatrix = importdata('input_data/pressure_adjoint_kWave_57600sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('input_data/pressure_adjoint_kWave_1600sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_adjoint_RT_1600sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_x_k.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_FISTA_tau5e18_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_SGD_27.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_FISTA_29.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_S-FISTA_10.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('S-FISTA, iter 10 - PROBLEM');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example72_S-FISTA_problem.fig');

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_PDHG_16.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);


%==============================
% RESULTS
%==============================
% Gradient
pixelPressureMatrix = importdata('results/pixelPressure_GD_tau5e18_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('GD - t = 5e18, l = 1e-2, iter = 50 - heterogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example72_GD_tau5e18_lambda1e-2_iter50.fig');

% Stochastic gradient
pixelPressureMatrix = importdata('results/pixelPressure_S-GD_tau5e18_lambda3e-4_epoch50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('S-GD - t = 5e18, l = 3e-4, iter = 50 - heterogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example72_S-GD_tau5e18_lambda1e-2_iter50.fig');

% FISTA
pixelPressureMatrix = importdata('results/pixelPressure_FISTA_tau5e18_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('FISTA - t = 5e18, l = 1e-2, iter = 50 - heterogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example72_FISTA_tau5e18_lambda1e-2_iter50.fig');

% Stochastic FISTA
%pixelPressureMatrix = importdata('results/pixelPressure_S-FISTA_tau1e18_lambda5e-4_epoch50.dat', ' ', 0);
%pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%plot_projection(pixelPressure, dx);
%a = axes;
%t = title('S-FISTA - t = 1e18, l = 1e-2, iter = 50 - heterogeneous SS');
%a.Visible = 'off'; 
%t.Visible = 'on'; 
%saveas(gcf, 'Example72_GD_tau1e18_lambda1e-2_iter50.fig');

% PDHG
pixelPressureMatrix = importdata('results/pixelPressure_PDHG_sigma1e0_tau5e18_theta1_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('PDHG - s = 1, t = 5e18, l = 1e-2, iter = 50 - heterogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example72_PDHG_sigma1_tau5e18_lambda1e-2_iter50.fig');


% Pixel pressure
pixelPressureMatrix = importdata('results/pixelPressure_adjoint_RT_57600sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
