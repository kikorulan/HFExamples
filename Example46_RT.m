%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex46_caustics;

close all;
%clear all;

load sensor_data.mat;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

% Make matrix symmetric
%c = ones(Nx, Ny);
%c(1:floor(Nx/2), :) = medium.sound_speed(1:floor(Nx/2), :);
%c(floor(Nx/2) + 1, :) = medium.sound_speed(floor(Nx/2), :);
%c(floor(Nx/2) + 2:end, :) = fliplr(medium.sound_speed(1:floor(Nx/2), :)')';
c = medium.sound_speed;
gridR.setCMatrix(c);


%========================================
% Impulse Response
%========================================
%%  % Set time
%%  dt = 2e-8;
%%  %dt = min(gridR.dx, gridR.dy)/c0/2;
%%  tMax = 2e-5;
%%  gridR.setTime(dt, tMax);
%%  % Compute impulse response
%%  gridR.impulse_additive('IV');
%%  save gridRT_impulse.mat gridR c0 dt tMax;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 800;% 800
nSources = 1;%256

% Parametrisation
tStep = dt;

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-1)/2*gridR.dx, 0);

% Sources
source(1) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
source(2) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
source(3) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
% Set initial pressure
gridR.setUMatrix(source_low.p0);
gridR.computeHamil(source(1));
gridR.setUMatrix(source_mid.p0);
gridR.computeHamil(source(2));
gridR.setUMatrix(source_high.p0);
gridR.computeHamil(source(3));

%%  clear source;
%%  source = gridR.computeForwardParallel(x, 0, pi, nRays, tStep, tMax, false);

% Save results
%save gridRT.mat gridR nRays nSources x -v7.3;
  
%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex46_caustics;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 375 600];
positionBar = [700 700 375 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Initial Pressure
%==============================
pressure = source_low.p0 + source_mid.p0 + source_high.p0;
figure;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
hold on;
surf(1e3*gridR.xAxis, 1e3*gridR.yAxis, pressure', 'EdgeColor', 'none');
view(2);
%title('Initial Pressure');
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'ytick', []);
%set(gca, 'yticklabel', []);
colorbar();
set(gcf, 'pos', positionYBar);
saveas(gcf, 'Example46_U', 'png'); 
saveas(gcf, 'Example46_U.fig'); 

%==============================
% Sound speed
%==============================
gridR.plot_soundSpeed();
set(gcf, 'pos', positionYBar);
colorbar();
saveas(gcf, 'Example46_C', 'png'); 
saveas(gcf, 'Example46_C.fig'); 

%==============================
% Traveled distance
%==============================
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      plot(gridR.tForward, real(source(1).xTD(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Traveled distance');
%%  xlabel('t (s)');

%==============================
% Sources and their respective Rays
%==============================
%%  for n = 1:nSources
%%     source(1).plot_rays(gridR, 10);
%%  end

%==============================
% Isotime curve
%==============================
%%  nCurves = 20;
%%  xIso = gridR.findIsoTime(source(1), 20, 1e-4, 1.2e-4);
%%  % Figure
%%  figure;
%%  flipGrid = 1./gridR.c';
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  hold on;
%%  surf(gridR.xAxis, gridR.yAxis, flipGrid, 'EdgeColor', 'none');
%%  hold on;
%%  for n = 1:nCurves
%%      plot3(xIso(n, :, 1), xIso(n, :, 2), repmat(1501, [1 size(xIso, 2)]), '-m');
%%  end
%%  for j = 1:20:nRays
%%      plot3(source(2).x(j, :, 1), source(2).x(j, :, 2), repmat(1501, [1 source(2).nPoints]), '-g');
%%  end
%%  view(2);
%%  legend('Isocurves for times tMin = 1e-4s, tMax = 1.2e-4s');

%==============================
% Sound speed + rays
%==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = nRays;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:nRays
%%      colourNum = floor(nColours*(j-1)/nRays) + 1;
%%      plot3(source(1).x(j, :, 1), source(1).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  %surf(gridR.xAxis, gridR.yAxis, 1./gridR.c', 'EdgeColor', 'none');
%%  imagesc(gridR.xAxis, gridR.yAxis, 1./gridR.c');
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', position);
%%  %saveas(gcf, 'Example46_rays_soundSpeed.fig');

%==============================
% Ray trajectories
%==============================
h = figure;
hold on;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
nColours = nRays;
colourList = summer(nColours);
[nRays nSteps dim] = size(source(1).x);
zVec = 1600*ones(1, nSteps);
for j = 1:5:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot3(1e3*source(1).x(j, :, 1), 1e3*source(1).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
end
xlabel('x [mm]');
ylabel('y [mm]');
box on;
set(gcf, 'pos', position);
saveas(gcf, 'Example46_rays', 'png');
saveas(gcf, 'Example46_rays.fig');


%==============================
% Amplitude
%==============================
%%  % Imagesc
%%  figure;
%%  surf(real(log(real(source(1).amplitude))), 'EdgeColor', 'none');
%%  view(2);
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('Amplitudes (log scale)');
%%  colorbar();
%%  hold on;
%%  %plot(1:1000, 350*ones(1, 1000), 'Color', 'r');
%%  %plot(1:1000, 450*ones(1, 1000), 'Color', 'r');
%%  %%  saveas(gcf, 'Example46_surf_amplitudeDin_log', 'png');
%%  %%  
%%  % Log plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      semilogy(gridR.tForward, real(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Amplitudes for rays with caustic (Gaussian Beam)');
%%  xlabel('t (s)');
%%  saveas(gcf, 'Example46_amplitudeGB_log', 'epsc');
%%  
%%  % Log plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      semilogy(gridR.tForward, imag(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Amplitudes for rays with caustic - imag');
%%  xlabel('t (s)');
%%  %saveas(gcf, 'Example46_amplitudeODE_log_caustic.fig');

%==============================
% Amplitude & Q - RT & GB
%==============================
%%  % Amplitude 
%%  figure;
%%  colours = winter(100);
%%  subplot(2, 1, 1);
%%  for n = 1:100
%%      plot(gridR.tForward, real(sourceRT(1).q(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  ax = gca;
%%  ax.GridAlpha = 0.5;
%%  grid on;
%%  box on;
%%  title('Q rays with caustic');
%%  xlabel('t (s)');
%%  subplot(2, 1, 2);
%%  for n = 1:100
%%      semilogy(gridR.tForward, real(sourceRT(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      axis([0 2e-5 1e-3 1]);
%%      hold on;
%%  end;
%%  ax = gca;
%%  ax.GridAlpha = 0.5;
%%  grid on;
%%  box on;
%%  title('Amplitude rays with caustic');
%%  xlabel('t (s)');
%%  set(gcf, 'pos', [700 700 500 500]);
%%  saveas(gcf, 'Example46_AQ', 'epsc');


% Amplitude Gaussian Beam
%%  figure;
%%  subplot(2, 1, 1);
%%  for n = 1:100
%%      plot(gridR.tForward, real(source(1).qGB(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  ax = gca;
%%  ax.GridAlpha = 0.5;
%%  grid on;
%%  box on;
%%  title('Q rays with caustic (Gaussian Beam)');
%%  xlabel('t (s)');
%%  subplot(2, 1, 2);
%%  for n = 1:100
%%      semilogy(gridR.tForward, real(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      axis([0 2e-5 1e-3 1]);
%%      hold on;
%%  end;
%%  ax = gca;
%%  ax.GridAlpha = 0.5;
%%  grid on;
%%  box on;
%%  title('Amplitudes rays with caustic (Gaussian Beam)');
%%  xlabel('t (s)');
%%  set(gcf, 'pos', [700 700 500 500]);
%%  saveas(gcf, 'Example46_AQ_GB', 'epsc');



%==============================
% Q
%==============================
%%  % Imagesc
%%  figure;
%%  imagesc(real(source(1).q));
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('Q');
%%  colorbar();
%%  %caxis([-8e-7 8e-7]);
%%  saveas(gcf, 'Example46_surf_qDin', 'png');

%hold on;
%plot(1:1000, 350*ones(1, 1000), 'Color', 'r');
%plot(1:1000, 450*ones(1, 1000), 'Color', 'r');
%saveas(gcf, 'Example46_qODE.fig');

%%  % Real plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      plot(gridR.tForward, real(source(1).q(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Qs for rays with caustic');
%%  xlabel('t (s)');
%%  saveas(gcf, 'Example46_qDin', 'epsc');
%%  
%%  % Imag plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      plot(gridR.tForward, imag(source(1).q(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Qs for rays with caustic - imag');
%%  xlabel('t (s)');
%%  %saveas(gcf, 'Example46_q.fig');

%==============================
% Q GB
%==============================
%%  % Imagesc
%%  figure;
%%  imagesc(real(source(1).qGB));
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('Q GB real');
%%  colorbar();
%%  
%%  
%%  figure;
%%  imagesc(imag(source(1).qGB));
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('Q GB imag');
%%  colorbar();



%caxis([-8e-7 8e-7]);


%==============================
% Keller-Maslov index
%==============================
%%  % Imagesc
%%  figure;
%%  imagesc(real(source(1).kIndex));
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('KM index real');
%%  colorbar();
%%  %caxis([-8e-7 8e-7]);
%%  
%%  % Imagesc
%%  figure;
%%  imagesc(imag(source(1).kIndex));
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('KM index imag');
%%  colorbar();
%%  %caxis([-8e-7 8e-7]);

%%  % Real plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      plot(gridR.tForward, real(source(1).Dx(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('K index for rays with caustic');
%%  xlabel('t (s)');

%%  %==============================
%%  % Curvature and Beamwidth
%%  %==============================
%%  
%%  % Curvature
%%  curvature = 1./source(1).n.*real(source(1).kIndex./source(1).q);
%%  figure;
%%  colours = cool(100);
%%  for n = 1:100
%%      plot(gridR.tForward, curvature(n+350, :), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Curvature');
%%  xlabel('t (s)');
%%  
%%  % Beamwidth
%%  beamwidth = sqrt(-1./imag(source(1).kIndex./source(1).q));
%%  figure;
%%  colours = cool(100);
%%  for n = 1:100
%%      plot(gridR.tForward, beamwidth(n+350, :), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('beamwidth');
%%  xlabel('t (s)');

%==============================
% Sound speed + ray selection
%==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = 320;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:nColours
%%      colourNum = j;
%%      plot3(source(1).x(j + 400-nColours/2, :, 1), source(1).x(j + 400-nColours/2, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  c = parula;
%%  c = flipud(c);
%%  colormap(c);
%%  imagesc(gridR.xAxis, gridR.yAxis, gridR.c');
%%  colorbar();
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', positionYBar);
%%  saveas(gcf, 'Example46_rays_soundSpeed_selection', 'png');
%%  saveas(gcf, 'Example46_rays_soundSpeed_selection.fig');

%==============================
% Time Signals
%==============================
%%  %source = sourceRT;
%%  source = sourceGB;
%%  % Norms
%%  normRT = max(real(source(1).aForward));
%%  normRT_2 = max(real(sourceRT(1).aForward));
%%  normKWave = max(sensor_data_low(1, :));
%%  figure;
%%  % Subplot 1
%%  subplot(2, 1, 1);
%%  hold on;
%%  axis([0 2e-5 -1.5 2]);
%%  grid on;
%%  box on;
%%  plot(gridR.tForward, source(1).aForward/normRT, 'Color', 'r', 'LineWidth', 2);
%%  plot(kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
%%  plot(gridR.tForward, source(2).aForward/normRT, 'Color', 'g', 'LineWidth', 2);
%%  plot(kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
%%  plot(gridR.tForward, source(3).aForward/normRT, 'Color', 'b', 'LineWidth', 2);
%%  plot(gridR.tForward, sourceRT(3).aForward/normRT_2, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
%%  plot(kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
%%  %legend('RT - bottom', 'kWave - bottom', 'RT - middle', 'kWave - middle', 'RT - top', 'kWave - top');
%%  legend('GB - bottom', 'kWave - bottom', 'GB - middle', 'kWave - middle', 'GB - top', 'RT - top', 'kWave - top');
%%  xlabel('t (s)');
%%  ylabel('Amplitude');
%%  title('Forward data - RT (Gaussian Beam) vs k-Wave');
%%  % Subplot 2
%%  subplot(2, 1, 2);
%%  hold on;
%%  axis([0 2e-5 -1 1]);
%%  grid on;
%%  box on;
%%  plot(gridR.tForward, source(1).aForward/normRT - sensor_data_low(1, :)/normKWave, 'Color', 'r', 'LineWidth', 2);
%%  plot(gridR.tForward, source(2).aForward/normRT - sensor_data_mid(1, :)/normKWave, 'Color', 'g', 'LineWidth', 2);
%%  plot(gridR.tForward, source(3).aForward/normRT - sensor_data_high(1, :)/normKWave, 'Color', 'b', 'LineWidth', 2);
%%  plot(gridR.tForward, sourceRT(3).aForward/normRT_2 - sensor_data_high(1, :)/normKWave, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
%%  legend('Error - bottom', 'Error - middle', 'Error - top', 'Error (RT) - top');
%%  xlabel('t (s)');
%%  ylabel('Amplitude');
%%  title('Error - RT ');
%%  set(gcf, 'pos', [700 700 1200 800]);
%%  saveas(gcf, 'Example46_signalsGB_error', 'png');
%%  saveas(gcf, 'Example46_signalsGB_error.fig');

%%  %==============================
%%  % Beam Signals
%%  %==============================
%%  % Norms
%%  normRT = max(real(source(1).aBeam));
%%  % Number of plots
%%  figure;
%%  hold on;
%%  grid on;
%%  box on;
%%  plot(gridR.tForward, real(source(1).aBeam)/normRT, 'Color', 'r', 'LineWidth', 2);
%%  plot(gridR.tForward, imag(source(1).aBeam)/normRT, 'Color', [0.8 0.2 0.2]);
%%  plot(gridR.tForward, real(source(2).aBeam)/normRT, 'Color', 'g', 'LineWidth', 2);
%%  plot(gridR.tForward, imag(source(2).aBeam)/normRT, 'Color', [0.2 0.8 0.2]);
%%  plot(gridR.tForward, real(source(3).aBeam)/normRT, 'Color', 'b', 'LineWidth', 2);
%%  plot(gridR.tForward, real(source(3).aBeam)/normRT, 'Color', [0.2 0.2 0.8]);
%%  legend('RT - real low', 'RT - imag low', 'RT - real mid', 'RT - imag mid', 'RT - real high', 'RT - imag high');
%%  xlabel('t (s)');
%%  ylabel('Amplitude');
%%  %title('Signals for focusing sound speed - ODE');

%%  %==============================
%%  % Time Signals
%%  %==============================
%%  % Norms
%%  normRT = max(imag(source(1).aForward));
%%  normKWave = max(sensor_data_low(1, :));
%%  % Number of plots
%%  figure;
%%  hold on;
%%  axis([0 2e-4 -1.5 2.5]);
%%  grid on;
%%  box on;
%%  plot(gridR.tForward, real(source(1).aForward)/normRT, 'Color', 'r', 'LineWidth', 2);
%%  plot(gridR.tForward, imag(source(1).aForward)/normRT, 'Color', 'r', 'LineWidth', 2);
%%  plot(kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
%%  plot(gridR.tForward, real(source(2).aForward)/normRT, 'Color', 'g', 'LineWidth', 2);
%%  plot(gridR.tForward, imag(source(2).aForward)/normRT, 'Color', 'g', 'LineWidth', 2);
%%  plot(kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
%%  plot(gridR.tForward, real(source(3).aForward)/normRT, 'Color', 'b', 'LineWidth', 2);
%%  plot(gridR.tForward, imag(source(3).aForward)/normRT, 'Color', 'b', 'LineWidth', 2);
%%  plot(kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
%%  legend('RT - low', 'kWave - low', 'RT - mid', 'kWave - mid', 'RT - high', 'kWave - high');
%%  xlabel('t (s)');
%%  ylabel('Amplitude');
%%  %title('Signals for focusing sound speed - ODE');
%%  %saveas(gcf, 'Example46_signalsODE_focus.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

