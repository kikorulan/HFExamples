cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex80_3D_veins_subsampled;
%cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex80_3D_veins_subsampled;

clear all;
close all;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('./input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);


%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data_kWave = h5read('./output_data/Example80_forward_output_25sensors.h5', '/p');
figure;
imagesc(sensor_data_kWave);

%==================================================
% TIME SIGNAL - RT
%==================================================
timeSignal = importdata(['./output_data/forwardSignal_25sensors_amplMod.dat'], ' ', 0);
sensor_data_RT = timeSignal(2:end, :);
figure;
imagesc(sensor_data_RT);

%==================================================
% COMPARE
%==================================================
% Normalize
normKW = max(sensor_data_kWave(:));
normRT = max(sensor_data_RT(:));
sensor_data_kWave_norm = sensor_data_kWave/normKW;
sensor_data_RT_norm = sensor_data_RT/normRT;

% Plot
index_sensor = 1;
figure;
plot(sensor_data_kWave_norm(index_sensor, :), 'Color', 'r');
hold on;
plot(sensor_data_RT_norm(index_sensor, :), 'Color', 'b');

%==================================================
% AMPLITUDE
%==================================================
amplitude     = importdata(['./output_data/Amplitude20.dat'], ' ', 0);
amplitude_mod = importdata(['./output_data/Amplitude20_amplMod.dat'], ' ', 0);
figure;
semilogy(amplitude(:, 1), 'Color', 'r')
hold on;
semilogy(amplitude_mod(:, 1), 'Color', 'b')
grid on;
ax = gca;
ax.GridAlpha = 0.4;


amplitude_norm = amplitude/amplitude(1);
amplitude_mod_norm = amplitude_mod/amplitude_mod(1);
figure;
semilogy(amplitude_norm(:, 1), 'Color', 'r', 'Linewidth', 1.5)
hold on;
semilogy(amplitude_mod_norm(:, 1), 'Color', 'b', 'Linewidth', 1.5)
grid on;
ax = gca;
ax.GridAlpha = 0.4;

