% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex77_3D_noisy_vessels;

clear all;
close all;


% Import data
filenameData = './input_data/forwardSignal_reference_14400sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

%==================================================
% ADD RANDOM NOISE
%==================================================
[nSensors, signalLength] = size(inputRT);
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
dlmwrite('./input_data/forwardSignal_reference_noisy5_14400sensors.dat', timeSignal_noisy, 'delimiter', ' ');
