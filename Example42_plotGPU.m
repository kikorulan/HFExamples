% Read data from files
cd ~/Documents/C++/HighFreq_3DRT/Build/bin/output_data/
%clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%==================================================
% Dimensions
%==================================================
% Import dimensions
filenameDimensions = 'dimensions.dat';
dim = importdata(filenameDimensions, delimiterIn, headerlinesIn);
Nx = dim(1, 1);
dx = dim(2, 1);
Ny = dim(1, 2);
dy = dim(2, 2);
Nz = dim(1, 3);
dz = dim(2, 3);

%%  %==================================================
%%  % GRADIENT
%%  %==================================================
%%  % Import data
%%  n_sound_speed = importdata('n.dat', delimiterIn, headerlinesIn);
%%  gradX = importdata('gradX.dat', delimiterIn, headerlinesIn);
%%  gradY = importdata('gradY.dat', delimiterIn, headerlinesIn);
%%  gradZ = importdata('gradZ.dat', delimiterIn, headerlinesIn);
%%  
%%  % Plot gradient
%%  figure;
%%  surf(n_sound_speed, 'EdgeColor', 'none');
%%  view(2); 
%%  title('N ');
%%  
%%  figure;
%%  surf(gradX, 'EdgeColor', 'none');
%%  view(2);
%%  title('Grad X');
%%  
%%  figure;
%%  surf(gradY, 'EdgeColor', 'none');
%%  view(2);
%%  title('Grad Y');
%%  
%%  figure;
%%  surf(gradZ, 'EdgeColor', 'none');
%%  view(2);
%%  title('Grad Z');

%%  %==================================================
%%  % TRAJECTORIES
%%  %==================================================
%%  % Import data
%%  filenameData = 'Trajectory0.dat';
%%  trajectories = importdata(filenameData, delimiterIn, headerlinesIn);
%%  
%%  % Read number of rays and steps
%%  [nSteps nRays] = size(trajectories);
%%  xCoord = trajectories(:, 1:3:nRays);
%%  yCoord = trajectories(:, 2:3:nRays);
%%  zCoord = trajectories(:, 3:3:nRays);
%%  
%%  nRays = nRays/3;
%%  % Plot the figure
%%  figure;
%%  ax = gca;
%%  ax.GridAlpha = 1;
%%  grid on;
%%  axis([0 Nx*dx 0 Ny*dy 0 Nz*dz]);
%%  hold on;
%%  colours = winter(nRays);
%%  for n = 1:nRays
%%      plot3(xCoord(:, n), yCoord(:, n), zCoord(:, n), 'Color', colours(n, :));
%%  end

%%  %==================================================
%%  % INITIAL PRESSURE
%%  %==================================================
%%  % Import data
%%  filenameData = 'InitialPressure0.dat';
%%  initialPressure = importdata(filenameData, delimiterIn, headerlinesIn);
%%  
%%  % Plot
%%  figure;
%%  surf(initialPressure, 'EdgeColor', 'none');
%%  view(2);

%%  %==================================================
%%  % AMPLITUDE
%%  %==================================================
%%  % Import data
%%  filenameData = 'Amplitude0.dat';
%%  amplitude = importdata(filenameData, delimiterIn, headerlinesIn);
%%  
%%  [nSteps nRays] = size(amplitude);
%%  % Plot
%%  figure;
%%  ax = gca;
%%  ax.GridAlpha = 1;
%%  grid on;
%%  colours = winter(nRays);
%%  for n = 1:nRays
%%      plot(amplitude(:, n), 'Color', colours(n, :));
%%      hold on;
%%  end

%==================================================
% TIME SIGNAL
%==================================================
% Import data
filenameData = 'ForwardSignal0.dat';
timeSignal = importdata(filenameData, delimiterIn, headerlinesIn);

% Plot
figure;
plot(timeSignal(:, 1), timeSignal(:, 2), 'Color', 'r');

%==================================================
% COMPARISON WITH K-WAVE
%==================================================
% Normalisation - RT
normRT = max(timeSignal(:, 2));
timeRT = timeSignal(:, 1);
signalRT = timeSignal(:, 2)/normRT;

% Normalisation - kWave
inputKWave = sensor_data_3balls_20p;
normKWave = max(inputKWave);
timeKWave = kgrid.t_array;
signalKWave = inputKWave/normKWave;

% Plot
figure;
grid on;
plot(timeRT, signalRT, 'Color', 'r', 'LineWidth', 2);
hold on;
plot(timeKWave, signalKWave, 'Color', 'blue', 'LineWidth', 2);
axis([0 1e-4 -1 1]);
legend('RT', 'k-Wave');
xlabel('time (s)');
ylabel('amplitude');
title('Comparison between RT and kWave');
saveas(gcf, 'Example42_RT_GPU.fig');
saveas(gcf, 'Example42_RT_GPU', 'png');

% Error
signalKWave_spline = spline(timeKWave, signalKWave, timeRT);
error = signalRT - signalKWave_spline;
figure;
grid on;
plot(timeRT, error, 'Color', 'r', 'LineWidth', 2);
legend('error');
xlabel('time (s)');
ylabel('error');
saveas(gcf, 'Example42_error_GPU.fig');
saveas(gcf, 'Example42_error_GPU', 'png');
