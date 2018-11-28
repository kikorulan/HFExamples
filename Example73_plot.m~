% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex70_3D_synchronization;

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

%============================================================================================================================================
% FORWARD PROBLEM
%============================================================================================================================================

%==================================================
% TIME SIGNAL - RT
%==================================================
% Import data
filenameData = 'output_data/ForwardSignal_delay_2_-1.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT_1 = timeSignal(2:end, :);

% Import data
filenameData = 'output_data/ForwardSignal_delay_1_0.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT_2 = timeSignal(2:end, :);

% Import data
filenameData = 'output_data/ForwardSignal_delay_1_0_delta.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT_3 = timeSignal(2:end, :);

%==================================================
% TIME SIGNAL - kWave
%==================================================
load ./input_data/sensor_data_veins_25sensors.mat;
sensor_data = h5read('output_data/Example70_forward_output_25sensors.h5', '/p');
inputKWave = sensor_data(:, :);

%==================================================
% Comparison
%==================================================
normRT_1 = max(inputRT_1(:));
normRT_2 = max(inputRT_2(:));
normRT_3 = max(inputRT_3(:));
normKW = max(inputKWave(:));

figure;
hold on;
plot(timeRT, inputRT_1(1, :)/normRT_1);
plot(timeRT, inputRT_2(1, :)/normRT_2);
plot(timeRT, inputRT_3(1, :)/normRT_3);
plot(kgrid.t_array, inputKWave(13, :)/normKW);
legend('RT 1', 'RT 2', 'RT 3', 'KW 1');

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================


%==================================================
% Reconstruction - kWave
%==================================================
p0_recon_PML = h5read('output_data/Example70_adjoint_output_1sensor.h5', '/p_final');
PML_size = 10;
p0_recon = p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
plot_projection(p0_recon, dx);
saveas(gcf, 'Example68_kWave_adjoint.fig');

%==================================================
% Reconstruction - RT
%==================================================
% Import data
pixelPressureMatrix = importdata('output_data/PixelPressure_delay_2_-1.dat', ' ', 0);
pixelPressure_1 = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure_1, dx);

% Import data
pixelPressureMatrix = importdata('output_data/PixelPressure_delay_1_0.dat', ' ', 0);
pixelPressure_2 = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure_2, dx);

% Import data
pixelPressureMatrix = importdata('output_data/PixelPressure_delay_1_0_delta.dat', ' ', 0);
pixelPressure_3 = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure_3, dx);


%==================================================
% Reconstruction - RT
%==================================================
maxKWave = max(p0_recon(:));
maxRT_1 = max(pixelPressure_1(:));
maxRT_2 = max(pixelPressure_2(:));
maxRT_3 = max(pixelPressure_3(:));
lineKWave = p0_recon(:, 121, 121);
lineRT_1 = pixelPressure_1(:, 121, 121);
lineRT_2 = pixelPressure_2(:, 121, 121);
lineRT_3 = pixelPressure_3(:, 121, 121);

figure;
hold on;
plot(lineRT_1/maxRT_1);
plot(lineRT_2/maxRT_2);
plot(lineRT_3/maxRT_3);
plot(lineKWave/maxKWave);
legend('RT 1', 'RT 2', 'KW');
