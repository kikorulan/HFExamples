%================================================================================
% Example: Study of caustics
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex48_RT;

close all;
%clear all;

load sensor_data.mat;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

% Make matrix symmetric
%c = ones(Nx, Ny);
%c(1:floor(Nx/2), :) = medium.sound_speed(1:floor(Nx/2), :);
%c(floor(Nx/2) + 1, :) = medium.sound_speed(floor(Nx/2), :);
%c(floor(Nx/2) + 2:end, :) = fliplr(medium.sound_speed(1:floor(Nx/2), :)')';
c = medium.sound_speed;
gridR.setCMatrix(c);

%%  %========================================
%%  % Impulse Response
%%  %========================================
%%  % Set time
%%  dt = 2e-7;
%%  %dt = min(gridR.dx, gridR.dy)/c0/2;
%%  tMax = 2.5e-4;
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
nRays = 100;% 800
nSources = 3;%256

% Parametrisation
tStep = dt;

% Sources locations
x = cat(3, floor((gridR.Nx-1)/2)*gridR.dx, 0);
x_low = cat(3, floor((gridR.Nx-1)/2)*gridR.dx, floor((gridR.Ny-1)/4)*gridR.dy);
x_mid = cat(3, floor((gridR.Nx-1)/2)*gridR.dx, floor((31*(gridR.Ny-1)/50))*gridR.dy);
x_high = cat(3, floor((gridR.Nx-1)/2)*gridR.dx, floor((4*(gridR.Ny-1)/5))*gridR.dy);

% Sources
source(1) = gridR.newSource(x, pi/2-1, pi/2+1, nRays, tStep, tMax);
source(2) = gridR.newSource(x, pi/2-1, pi/2+1, nRays, tStep, tMax);
source(3) = gridR.newSource(x, pi/2-1, pi/2+1, nRays, tStep, tMax);
% Set initial pressure

% Source 0
gridR.setUMatrix(source_low.p0);
gridR.forward_trajectories(source(1));
gridR.forward_amplitudeGB(source(1));
gridR.forward_revAmplitude(source(1));
gridR.forward_pressure(source(1));
gridR.forward_subrays(source(1), 10, gridR.dx);
gridR.forward_beam(source(1));
gridR.forward_timeSignal(source(1));

gridR.setUMatrix(source_mid.p0);
gridR.forward_trajectories(source(2));
gridR.forward_amplitudeGB(source(2));
gridR.forward_revAmplitude(source(2));
gridR.forward_pressure(source(2));
gridR.forward_subrays(source(2), 10, gridR.dx);
gridR.forward_beam(source(2));
gridR.forward_timeSignal(source(2));

gridR.setUMatrix(source_high.p0);
gridR.forward_trajectories(source(3));
gridR.forward_amplitudeGB(source(3));
gridR.forward_revAmplitude(source(3));
gridR.forward_pressure(source(3));
gridR.forward_subrays(source(3), 10, gridR.dx);
gridR.forward_beam(source(3));
gridR.forward_timeSignal(source(3));


%%  % Subrays
%%  gridR.forward_subrays(source(2), 5, gridR.dx);
%%  gridR.forward_subrays(source(3), 5, gridR.dx);

%%  % Pixel time
%%  gridR.forward_raysToGrid_multiray(source(1));
%%  gridR.forward_raysToGrid_multiray(source(2));
%%  gridR.forward_raysToGrid_multiray(source(3));

%%  clear source;
%%  source = gridR.computeForwardParallel(x, 0, pi, nRays, tStep, tMax, false);

% Save results
%save gridRT.mat gridR nRays nSources x -v7.3;
  
%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex48_RT;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%%  %==============================
%%  % Initial Pressure
%%  %==============================
%%  pressure = source_low.p0 + source_mid.p0 + source_high.p0;
%%  figure;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  hold on;
%%  surf(gridR.xAxis, gridR.yAxis, pressure', 'EdgeColor', 'none');
%%  view(2);
%%  title('Initial Pressure');
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  set(gcf, 'pos', position);
%%  %saveas(gcf, 'Example46_U.fig'); 

%==============================
% Sound speed
%==============================
gridR.plot_soundSpeed();
set(gcf, 'pos', positionYBar);
colorbar();
saveas(gcf, 'Example48_C', 'epsc'); 

%%  %==============================
%%  % Traveled distance
%%  %==============================
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
%%      source(n).plot_rays(gridR, 100);
%%  end

h = source(1).plot_subrays(gridR, 50);
saveas(gcf, 'Example48_subrays', 'epsc');

%%  %==============================
%%  % Sound speed + rays
%%  %==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = nRays;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:nRays
%%      colourNum = floor(nColours*(j-1)/nRays) + 1;
%%      plot3(source(2).x(j, :, 1), source(2).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  %surf(gridR.xAxis, gridR.yAxis, 1./gridR.c', 'EdgeColor', 'none');
%%  imagesc(gridR.xAxis, gridR.yAxis, 1./gridR.c');
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', position);
%%  %saveas(gcf, 'Example46_rays_soundSpeed.fig');

%==============================
% Amplitude
%==============================
source(1).plot_amplitude(1, true);
h1 = source(1).plot_amplitude_subrays(1, true);
saveas(h1, 'Example48_amplitude_subrays', 'epsc');

%==============================
% Q
%==============================
source(1).plot_q(1, true);

%==============================
% Pressure
%==============================
source(1).plot_pressure();
[~, h] = source(1).plot_pressure_subrays(50);

%==============================
% Keller-Maslov index
%==============================

% Real plot
figure;
nRays = size(source(1).kIndex, 1);
colours = winter(nRays);
for n = 1:nRays
    plot(gridR.tForward, real(source(1).kIndex(n, :)), 'Color', colours(n, :));
    hold on;
end;
grid on;
box on;
title('K index for rays with caustic');
xlabel('t (s)');

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


%%  %==============================
%%  % Sound speed + ray selection
%%  %==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = 100;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:100
%%      colourNum = j;
%%      plot3(source(1).x(j + 350, :, 1), source(1).x(j + 350, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  surf(gridR.xAxis, gridR.yAxis, 1./gridR.c', 'EdgeColor', 'none');
%%  %imagesc(gridR.xAxis, gridR.yAxis, 1./gridR.c');
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', position);
%%  %saveas(gcf, 'Example46_rays_soundSpeed_selection.fig');



%==============================
% Time Signals
%==============================
% Norms
normRT = max(real(source(1).aForward));
normKWave = max(sensor_data_low(1, :));
% Number of plots
figure;
hold on;
axis([0 2e-4 -1.5 2.5]);
grid on;
box on;
plot(gridR.tForward, real(source(1).aForward)/normRT, 'Color', 'r', 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_low(1, :)/normKWave, 'Color', [0.8 0.2 0.2]);
plot(gridR.tForward, real(source(2).aForward)/normRT, 'Color', 'g', 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_mid(1, :)/normKWave, 'Color', [0.2 0.8 0.2]);
plot(gridR.tForward, source(3).aForward/normRT, 'Color', 'b', 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_high(1, :)/normKWave, 'Color', [0.2 0.2 0.8]);
legend('RT - low', 'kWave - low', 'RT - mid', 'kWave - mid', 'RT - high', 'kWave - high');
xlabel('t (s)');
ylabel('Amplitude');
%title('Signals for focusing sound speed - Gaussian Beam');
saveas(gcf, 'Example48_signalsGB_focus', 'epsc');

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

%%  %==============================
%%  % Time Signals
%%  %==============================
%%  pixelAttenuation_low = permute(fliplr(real(source(1).pixelAttenuation)), [2 1 3]);
%%  scrollView(pixelAttenuation_low, 3, [-1 1]);

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

