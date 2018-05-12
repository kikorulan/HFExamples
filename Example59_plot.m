% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex59_RT_3D;
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
%==============================
% SOUND SPEED
%==============================
%%  sound_speed_matrix = importdata('input_data/sound_speed.dat', delimiterIn, headerlinesIn);
%%  sound_speed = matrix2cube(sound_speed_matrix, Nz);
%%  plot_cube(sound_speed);

%==============================
% FILTERS
%==============================
%%  filters = importdata('output_data/Filters.dat', delimiterIn, headerlinesIn);
%%  nColors = cool(100);
%%  figure;
%%  hold on;
%%  for i = 1:100
%%      plot(filters(i, :), 'Color', nColors(i, :));
%%  end

%==================================================
% INITIAL PRESSURE
%==================================================
% Import data
filenameData = 'input_data/initial_pressure_veins.dat';
initial_pressure_matrix = importdata(filenameData, delimiterIn, headerlinesIn);
initial_pressure = matrix2cube(initial_pressure_matrix, Nz);

[h, h1, h2, hlink] = plot_pixel(inputKWave_norm, 1, 1, false, initial_pressure);

%==================================================
% AMPLITUDE
%==================================================
%%  for i = 0:25
%%      % Import data
%%      filenameData = ['output_data/Amplitude' int2str(i) '.dat'];
%%      amplitude = importdata(filenameData, delimiterIn, headerlinesIn);
%%      % Plot
%%      figure;
%%      surf(abs(amplitude), 'EdgeColor', 'none');
%%      view(2);
%%      pause
%%  end

%==================================================
% TIME SIGNAL - RT
%==================================================
% Import data
filenameData = 'output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, delimiterIn, headerlinesIn);
% Plot
figure;
imagesc(timeSignal(2:end, :));
box on;

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('output_data/Example59_forward_output.h5', '/p');
% Plot
figure;
imagesc(sensor_data);
box on;


%==================================================
% COMPARISON WITH K-WAVE
%==================================================
% Import data
load input_data/sensor_data_veins;
sensor_data = h5read('output_data/Example59_forward_output.h5', '/p');
filenameData = 'output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, ' ', 0);
% Input
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
timeKWave = kgrid.t_array;
inputKWave = sensor_data;
nSensors = size(sensor_data, 1);


normRT = max(inputRT(:));
normKWave = max(inputKWave(:));

% Norm difference
difNorm = norm(inputRT(:)/normRT - inputKWave(:)/normKWave);

% Plot all sensors
for i = 1:100:nSensors
    % Normalisation - RT
    %normRT = max(inputRT(i, :));
    signalRT = inputRT(i, :)/normRT;

    % Normalisation - kWave
    %normKWave = max(inputKWave(i, :));
    signalKWave = inputKWave(i, :)/normKWave;
    
    % Plot
    figure;
    grid on;
    plot(timeRT, signalRT, 'Color', 'r', 'LineWidth', 2);
    hold on;
    plot(timeKWave, signalKWave, 'Color', 'blue', 'LineWidth', 2);
    axis([0 1.5e-5 -1 1]);
    legend('RT', 'k-Wave');
    xlabel('time (s)');
    ylabel('amplitude');
    grid on;
    title(['RT vs kWave - ', int2str(i)]);
    %saveas(gcf, ['output_data/Example59_forwardSignal_sensor', int2str(i), '.fig']);
    %saveas(gcf, ['output_data/Example59_forwardSignal_sensor', int2str(i)], 'png');
end

%%  for i = 1:10:nSensors
%%      % Error
%%      signalRT = inputRT(i, :)/normRT;
%%      signalKWave = inputKWave(i, :)/normKWave;
%%      % Spline
%%      signalKWave_spline = spline(timeKWave, signalKWave, timeRT);
%%      error = signalRT - signalKWave_spline;
%%      figure;
%%      grid on;
%%      plot(timeRT, error, 'Color', 'r', 'LineWidth', 2);
%%      axis([0 1.5e-5 -.2 .2]);
%%      legend('error');
%%      xlabel('time (s)');
%%      ylabel('error');
%%      grid on;
%%      title(['RT vs kWave - ', int2str(i)]);
%%      %saveas(gcf, ['output_data/Example59_forwardSignal_error', int2str(i), '.fig']);
%%      %saveas(gcf, ['output_data/Example59_forwardSignal_error', int2str(i)], 'png');
%%  end
%==========================================================================================
% INVERSE PROBLEM
%==========================================================================================
%%  % Position
%%  position     = [700 700 300 630];
%%  positionY    = [700 700 320 630];
%%  positionBar  = [700 700 363 630];
%%  positionYBar = [700 700 390 630];

%==================================================
% IMPORT DATA
%==================================================
%%  % Attenuation
%%  pixelAttenuationMatrix = importdata('output_data/PixelAttenuation0.dat', delimiterIn, headerlinesIn);
%%  pixelAttenuation = matrix2cube(pixelAttenuationMatrix, Nz);
%%  % Propagation time
%%  pixelTimeMatrix = importdata('output_data/PixelTime0.dat', delimiterIn, headerlinesIn);
%%  pixelTime = matrix2cube(pixelTimeMatrix, Nz);
%%  % Reconstruction by sensor
%%  for i = 1:26
%%      % RT
%%      p0_matrix = importdata(['output_data/PixelPressure_' int2str(i-1) '.dat'], delimiterIn, headerlinesIn);
%%      p0_recon_RT{i} = matrix2cube(p0_matrix, Nz);
%%      p0_recon_RT{i}(p0_recon_RT{i} == inf) = 0;
%%      p0_recon_RT{i}(p0_recon_RT{i} == -inf) = 0;
%%      % kWave
%%      p0_recon_PML = h5read(['output_data/Example56_adjoint_output_' int2str(i) '.h5'], '/p_final');
%%      PML_size = 10;
%%      p0_recon_KW{i} = p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
%%  end

%==================================================
% ATTENUATION
%==================================================
%%  % Plot
%%  figure;
%%  title('Attenuation');
%%  surf(pixelAttenuation(:, :, 64), 'EdgeColor', 'none');
%%  view(2);

%==================================================
% PROPAGATION TIME
%==================================================
%%  % Plot
%%  figure;
%%  title('Propagation time');
%%  surf(pixelTime(:, :, 64), 'EdgeColor', 'none');
%%  view(2);

%==================================================
% Reconstruction by sensor - RT/kWave
%==================================================
%%  pixelPressureMatrix = zeros(16384, 128);
%%  for i = 1:26
%%      % Choose slice
%%      inputRT = p0_recon_RT{i}(:, :, 64);
%%      normRT = max(inputRT(:));
%%      % Plot
%%      h1 = figure;
%%      surf(x_axis, y_axis, inputRT, 'EdgeColor', 'none');
%%      axis([0 x_axis(end) 0 y_axis(end)]);
%%      view(2);
%%      colorbar();
%%      xlabel('y (m)');
%%      ylabel('x (m)');
%%      title(['RT recon - z = 0.063, sensor = ' int2str(i)]);
%%  
%%      % Choose slice
%%      inputKWave = p0_recon_KW{i}(:, :, 64);    
%%      normKWave = max(inputKWave(:));
%%      % Plot figure
%%      h2 = figure;
%%      surf(x_axis, y_axis, inputKWave , 'EdgeColor', 'none');
%%      axis([0 x_axis(end) 0 y_axis(end)]);
%%      view(2);
%%      colorbar();
%%      xlabel('y (m)');
%%      ylabel('x (m)');
%%      title(['kWave recon - z = 0.063, sensor = ' int2str(i+1)]);
%%  
%%      % Error
%%      errorRecon = inputRT/normRT - inputKWave/normKWave;
%%      h3 = figure;
%%      surf(x_axis, y_axis, errorRecon, 'EdgeColor', 'none');
%%      axis([0 x_axis(end) 0 y_axis(end)]);
%%      view(2);
%%      colorbar();
%%  
%%      pause;
%%      close(h1); close(h2); close(h3);
%%  
%%  end

%==================================================
% Reconstruction - RT
%==================================================
% Import data
pixelPressureMatrix = importdata('output_data/PixelPressure.dat', delimiterIn, headerlinesIn);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_pixel(initial_pressure, 10, 1, pixelPressure);

%==================================================
% Reconstruction - kWave
%==================================================
p0_recon_PML = h5read('output_data/Example59_adjoint_output.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
plot_pixel(p0_recon, 10, 1);

%==================================================
% Comparison
%==================================================
%%  reconDif = inputRT_norm - p0_recon/normKW;
%%  % Plot figure
%%  scrollView(permute(fliplr(permute(real(reconDif), [2 1 3])), [2 1 3]), 3, [-0.2 .2])
%%  xlabel('y (m)');
%%  ylabel('x (m)');
%%  title('kWave recon - z = 0.063');
%%  saveas(gcf, 'output_data/Example57_error_recon.fig');


%==================================================
% COMPARISON
%==================================================
%%  % Plot
%%  figure;
%%  surf(x_axis, y_axis, inputRT_norm-inputKWave_norm, 'EdgeColor', 'none');
%%  axis([0 x_axis(end) 0 y_axis(end)]);
%%  view(2);
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Error RT - kWave - z = 0.063');
%%  saveas(gcf, 'Example56_error_recon', 'png');
%%  saveas(gcf, 'Example56_error_recon.fig');

%==================================================
% ENERGY ANALYSIS
%==================================================
%%  maxRT = 0;
%%  maxKW = 0;
%%  nSensors = 26;
%%  for i = 1:nSensors
%%      energy_RT{i} = permute(sqrt(sum(sum((p0_recon_RT{i}).^2, 1), 2)), [3 2 1]);
%%      energy_KW{i} = permute(sqrt(sum(sum((p0_recon_KW{i}).^2, 1), 2)), [3 2 1]);
%%      maxRT = max([energy_RT{i}; maxRT]);
%%      maxKW = max([energy_KW{i}; maxKW]);
%%  end
%%  normRT = max(energy_RT{4});
%%  normKW = max(energy_KW{4});
%%  for i = 1:nSensors
%%      energy_RT{i} = energy_RT{i}/normRT;
%%      energy_KW{i} = energy_KW{i}/normKW;
%%  end
%%  
%%  % kWave
%%  figure;
%%  nColor = cool(nSensors);
%%  for i = 1:nSensors
%%      hold on;
%%      plot(energy_KW{i}, 'Color', nColor(i, :));
%%  end
%%  
%%  % RT
%%  figure;
%%  nColor = cool(nSensors);
%%  for i = 1:nSensors
%%      hold on;
%%      plot(energy_RT{i}, 'Color', nColor(i, :));
%%  end
%%  
%%  % Sensor by sensor
%%  for i = 1:nSensors
%%      figure;
%%      plot(energy_RT{i}, 'Color', 'g', 'LineWidth', 2);
%%      hold on;
%%      plot(energy_KW{i}, 'Color', 'b', 'LineWidth', 2);
%%      %plot(energy_RT{i} - energy_KW{i}, 'Color', 'b', 'LineWidth', 2);
%%      title(['Energy for Sensor ' num2str(i)]);
%%      %legend('Ray Tracing', 'kWave', 'Difference');
%%  end

%==================================================
% COLUMN ANALYSIS
%==================================================
%%  % Colormap
%%  run colourMap;
%%  nSensors = 26;
%%  % Points
%%  p1x = round(3*Nx/4);
%%  p1y = round(1*Ny/4);
%%  p2x = round(1*Nx/4);
%%  p2y = round(2*Ny/4);
%%  p3x = round(3*Nx/4);
%%  p3y = round(3*Ny/4);
%%  % Select signals
%%  for i = 1:nSensors
%%      column1_RT{i} = permute(p0_recon_RT{i}(p1x, p1y, :), [3 2 1]);
%%      column2_RT{i} = permute(p0_recon_RT{i}(p2x, p2y, :), [3 2 1]);
%%      column3_RT{i} = permute(p0_recon_RT{i}(p3x, p3y, :), [3 2 1]);
%%      column1_KW{i} = permute(p0_recon_KW{i}(p1x, p1y, :), [3 2 1]);
%%      column2_KW{i} = permute(p0_recon_KW{i}(p2x, p2y, :), [3 2 1]);
%%      column3_KW{i} = permute(p0_recon_KW{i}(p3x, p3y, :), [3 2 1]);
%%  end
%%  
%%  % Plot
%%  for i = 1:nSensors
%%      % Normalize signals
%%      normColRT = max(column2_RT{i});
%%      normColKW = max(column2_KW{i});
%%  
%%      % Plot
%%      figure;
%%      plot(column1_RT{i}/normColRT, 'Color', colourMapV(1), 'LineWidth', 2); hold on;
%%      plot(column2_RT{i}/normColRT, 'Color', colourMapV(2), 'LineWidth', 2);
%%      plot(column3_RT{i}/normColRT, 'Color', colourMapV(3), 'LineWidth', 2);
%%      plot(column1_KW{i}/normColKW, 'Color', colourMapV(4));
%%      plot(column2_KW{i}/normColKW, 'Color', colourMapV(5));
%%      plot(column3_KW{i}/normColKW, 'Color', colourMapV(6));
%%      legend('RT 1', 'RT 2', 'RT 3', 'KW 1', 'KW 2', 'KW 3');
%%      title(['Sensor ' int2str(i)]);
%%  end


%%  %==================================================
%%  % VIDEO
%%  %==================================================
%%  % Import data
%%  pixelPressureMatrix = importdata('output_data/PixelPressure.dat', delimiterIn, headerlinesIn);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  p0_recon_PML = h5read('output_data/Example56_adjoint_output.h5', '/p_final');
%%  PML_size = 10;
%%  p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
%%  % Normalisation - RT
%%  inputRT = pixelPressure;
%%  inputRT(isinf(inputRT)) = 0;
%%  sortedNormRT = sort(inputRT(:));
%%  normRT = sortedNormRT(end-floor(length(sortedNormRT)/1e4));
%%  inputRT_norm = inputRT/normRT;
%%  % Normalisation - kWave
%%  inputKWave = p0_recon;
%%  sortedNormKW = sort(inputKWave(:));
%%  normKWave = sortedNormKW(end-floor(length(sortedNormKW)/1e4));
%%  inputKWave_norm = inputKWave/normKWave;
%%  
%%  % Plot
%%  fig = figure;
%%  fig.Units = 'points';
%%  fig.Position(3:4) = [1000 500];
%%  fig.PaperUnits = 'points';
%%  fig.PaperPosition(3:4) = [1000 500];
%%  ax11 = axes('units', 'points', 'position', [  50, 50, 400, 400]);
%%  ax12 = axes('units', 'points', 'position', [ 500, 50, 450, 400]);
%%  for i = 1:Nz
%%      subplot(ax11);
%%      surf(x_axis, y_axis, inputRT_norm(:, :, i), 'EdgeColor', 'none');
%%      axis([0 x_axis(end) 0 y_axis(end)]);
%%      view(2);
%%      caxis([0 1]);
%%      title(['RT reconstruction - z = ' num2str(i)])
%%      xlabel('y [m]');
%%      ylabel('x [m]');
%%      subplot(ax12);
%%      surf(x_axis, y_axis, inputKWave_norm(:, :, i), 'EdgeColor', 'none');
%%      axis([0 x_axis(end) 0 y_axis(end)]);
%%      view(2);
%%      caxis([0 1]);
%%      title(['kWave reconstruction - z = ' num2str(i)])
%%      xlabel('y [m]');
%%      ylabel('x [m]');
%%      colorbar();
%%      saveas(gcf, ['Example56_recon_' num2str(i)], 'jpg');
%%  end
%%  
%%  
%%  % Create New Video with the Image Sequence
%%  outputVideo = VideoWriter('Reconstruction.avi');
%%  outputVideo.FrameRate = 5;
%%  open(outputVideo)
%%  
%%  % Loop through the image sequence, load each image, and then write it to the video.
%%   for i = 1:Nz
%%     img = imread(['Example56_recon_' num2str(i) '.jpg']);
%%     writeVideo(outputVideo,img)
%%  end
%%  close(outputVideo)



%==================================================
% SIGNALS after convolution
%==================================================
%%  signalConvolve = importdata('SignalConvolve0.dat', delimiterIn, headerlinesIn);
%%  [nFilters lFilter] = size(signalConvolve);
%%  
%%  figure;
%%  hold on; 
%%  colours = winter(nFilters);
%%  for n = 1:nFilters
%%      plot(signalConvolve(n, :), 'Color', colours(n, :));    
%%  end

%==================================================
% FILTERS
%==================================================
%%  filters = importdata('Filters.dat', delimiterIn, headerlinesIn);
%%  [nFilters lFilter] = size(filters);
%%  
%%  figure;
%%  hold on;
%%  colours = winter(nFilters);
%%  for n = 1:nFilters
%%      plot(filters(n, :), 'Color', colours(n, :));    
%%  end

