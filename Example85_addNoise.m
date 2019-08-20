% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex85_3D_veins_subsampled;
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled;

clear all;
close all;


% Import data
filenameData = './input_data/forwardSignal_reference_3600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

% Dimensions
[nSensors, signalLength] = size(inputRT);

%==================================================
% ADD RANDOM NOISE - Multiplicative
%==================================================
% Random number generator
rng(1);
noise = normrnd(0,0.05,[nSensors,signalLength]);
% Add noise
inputRT_noisy = inputRT + inputRT.*noise;

% Plot noise
figure;
plot(noise(1, :));
% Plot input
figure;
plot(inputRT(1, :), 'LineWidth', 1.5);
hold on;
plot(inputRT_noisy(1, :), 'Color', 'r');
% Write output
timeSignal_noisy = [timeRT; inputRT_noisy];
dlmwrite('./input_data/forwardSignal_reference_noisy5_3600sensors.dat', timeSignal_noisy, 'delimiter', ' ');

%==================================================
% ADD RANDOM NOISE - Additive
%==================================================
% Random number generator
rng(1);
noise = normrnd(0,5e-3,[nSensors,signalLength]);
% Add noise
inputRT_noisy = inputRT + max(inputRT(:)).*noise;

% Plot noise
figure;
plot(noise(1, :));
% Plot input
figure;
plot(inputRT(1, :), 'LineWidth', 1.5);
hold on;
plot(inputRT_noisy(1, :), 'Color', 'r');

% Plot input and input noisy
figure; imagesc(inputRT(1:60, :));
figure; imagesc(inputRT_noisy(1:60, :));
% Write output
timeSignal_noisy = [timeRT; inputRT_noisy];
dlmwrite('./input_data/forwardSignal_reference_noisyAdd5_3600sensors.dat', timeSignal_noisy, 'delimiter', ' ');
