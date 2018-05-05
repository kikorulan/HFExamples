% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex52_RT_3D
%clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%==================================================
% Dimensions
%==================================================
% Import dimensions
filenameDimensions = 'input_data/dimensions.dat';
dim = importdata(filenameDimensions, delimiterIn, headerlinesIn);
Nx = dim(1, 1);
dx = dim(2, 1);
Ny = dim(1, 2);
dy = dim(2, 2);
Nz = dim(1, 3);
dz = dim(2, 3);

% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;

%==========================================================================================
% FORWARD PROBLEM
%==========================================================================================

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
%%  % Plot
%%  figure;
%%  surf(log(abs(amplitude)), 'EdgeColor', 'none');
%%  view(2);

%%  %==================================================
%%  % TIME SIGNAL
%%  %==================================================
%%  % Import data
%%  filenameData = 'output_data/ForwardSignal0.dat';
%%  timeSignal = importdata(filenameData, delimiterIn, headerlinesIn);
%%  
%%  % Plot
%%  figure;
%%  plot(timeSignal(:, 1), timeSignal(:, 2), 'Color', 'r');
%%  
%%  %==================================================
%%  % COMPARISON WITH K-WAVE
%%  %==================================================
%%  % Normalisation - RT
%%  normRT = max(timeSignal(:, 2));
%%  timeRT = timeSignal(:, 1);
%%  signalRT = timeSignal(:, 2)/normRT;
%%  
%%  % Normalisation - kWave
%%  inputKWave = sensor_data_3balls;
%%  normKWave = max(inputKWave);
%%  timeKWave = kgrid.t_array;
%%  signalKWave = inputKWave/normKWave;
%%  
%%  % Plot
%%  figure;
%%  grid on;
%%  plot(timeRT, signalRT, 'Color', 'r', 'LineWidth', 2);
%%  hold on;
%%  plot(timeKWave, signalKWave, 'Color', 'blue', 'LineWidth', 2);
%%  axis([0 1e-4 -1 1]);
%%  legend('RT', 'k-Wave');
%%  xlabel('time (s)');
%%  ylabel('amplitude');
%%  title('Comparison between RT and kWave');
%%  saveas(gcf, 'Example52_RT.fig');
%%  saveas(gcf, 'Example52_RT', 'png');
%%  
%%  % Error
%%  signalKWave_spline = spline(timeKWave, signalKWave, timeRT);
%%  error = signalRT - signalKWave_spline;
%%  figure;
%%  grid on;
%%  plot(timeRT, error, 'Color', 'r', 'LineWidth', 2);
%%  legend('error');
%%  xlabel('time (s)');
%%  ylabel('error');
%%  saveas(gcf, 'Example52_error.fig');
%%  saveas(gcf, 'Example52_error', 'png');

%==========================================================================================
% INVERSE PROBLEM
%==========================================================================================

%%  %==================================================
%%  % ATTENUATION
%%  %==================================================
%%  % Import data
%%  pixelAttenuationMatrix = importdata('output_data/PixelAttenuation0.dat', delimiterIn, headerlinesIn);
%%  pixelAttenuation = matrix2cube(pixelAttenuationMatrix, Nz);
%%  % Plot
%%  figure;
%%  title('Attenuation');
%%  surf(pixelAttenuation(:, :, 64), 'EdgeColor', 'none');
%%  view(2);
%%  
%%  %==================================================
%%  % PROPAGATION TIME
%%  %==================================================
%%  % Import data
%%  pixelTimeMatrix = importdata('output_data/PixelTime0.dat', delimiterIn, headerlinesIn);
%%  pixelTime = matrix2cube(pixelTimeMatrix, Nz);
%%  % Plot
%%  figure;
%%  title('Propagation time');
%%  surf(pixelTime(:, :, 64), 'EdgeColor', 'none');
%%  view(2);

%==================================================
% Reconstruction - RT
%==================================================
% Import data
pixelPressureMatrix = importdata('output_data/PixelPressure.dat', delimiterIn, headerlinesIn);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
% Normalisation - RT
inputRT = pixelPressure(:, :, 64);
inputRT(isinf(inputRT)) = 0;
normRT = max(inputRT(:));
inputRT_norm = inputRT/normRT;
% Plot
figure;
surf(x_axis, y_axis, inputRT_norm, 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('RT recon - z = 0.063');
saveas(gcf, 'Example52_RT_recon', 'png');
saveas(gcf, 'Example52_RT_recon.fig');

%==============================
% Reconstruction - kWave
%==============================
p0_recon_PML = h5read('output_data/Example52_adjoint_output.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));

% Normalisation - kWave
inputKWave = p0_recon(:, :, 64);
normKWave = max(inputKWave(:));
inputKWave_norm = inputKWave/normKWave;

% Plot figure
figure;
surf(x_axis, y_axis, inputKWave_norm, 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('kWave recon - z = 0.063');
saveas(gcf, 'Example52_kWave_recon', 'png');
saveas(gcf, 'Example52_kWave_recon.fig');

%==================================================
% COMPARISON
%==================================================
% Plot
figure;
surf(x_axis, y_axis, inputRT_norm-inputKWave_norm, 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
xlabel('x (m)');
ylabel('y (m)');
title('Error RT - kWave - z = 0.063');
saveas(gcf, 'Example52_error_recon', 'png');
saveas(gcf, 'Example52_error_recon.fig');

%%  %==================================================
%%  % SIGNALS after convolution
%%  %==================================================
%%  signalConvolve = importdata('SignalConvolve0.dat', delimiterIn, headerlinesIn);
%%  [nFilters lFilter] = size(signalConvolve);
%%  
%%  figure;
%%  hold on;
%%  colours = winter(nFilters);
%%  for n = 1:nFilters
%%      plot(signalConvolve(n, :), 'Color', colours(n, :));    
%%  end

%%  %==================================================
%%  % FILTERS
%%  %==================================================
%%  filters = importdata('Filters.dat', delimiterIn, headerlinesIn);
%%  [nFilters lFilter] = size(filters);
%%  
%%  figure;
%%  hold on;
%%  colours = winter(nFilters);
%%  for n = 1:nFilters
%%      plot(filters(n, :), 'Color', colours(n, :));    
%%  end

