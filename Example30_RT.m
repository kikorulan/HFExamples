%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex30_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex30_interf;

close all;
%clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
load sensor_data.mat;
%========================================
% Grid definition
%========================================
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]

Rgrid = gridRT(Nx, dx, Ny, dy);
% Time signal
dt = 3e-7;
tMax = 4e-4;
Rgrid.setTime(dt, tMax);

% Sound speed
c0 = 1500;
Rgrid.setCMatrix(medium.sound_speed);

% Build initial pressure
Rgrid.setUMatrix(sourceKW.p0);

% Impulse Response
Rgrid.impulse_additive('IV');

%========================================
% Ray Shooting
%========================================

% Number of rays & sources
nRays = 1000;
nSources = 3;

tMax = sqrt((Nx*dx)^2+(Ny*dy)^2)/c0;
tStep = dt;

% Sources locations
x1(1, 1, 1) = 2*Nx*dx/3;
x1(1, 1, 2) = 0; % Source 1
x2(1, 1, 1) = Nx*dx/4;
x2(1, 1, 2) = (Ny-1)*dy; % Source 2
x3(1, 1, 1) = 3*Nx*dx/4; 
x3(1, 1, 2) = (Ny-1)*dy; % Source 3

% Create the new sources
Rgrid.setDeltaX(Rgrid.dx/100);
source(1) = Rgrid.newSource(x1, 0, pi, nRays, tStep, tMax);
source(2) = Rgrid.newSource(x2, 0, -pi, nRays, tStep, tMax);
source(3) = Rgrid.newSource(x3, 0, -pi, nRays, tStep, tMax);

% Sources
for n = 1:nSources
    Rgrid.computeHamil(source(n));
end

% Save results
save gridRT.mat Rgrid source nRays nSources x1 x2 x3 -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex30_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex30_interf;

load gridRT.mat;
load sensor_data.mat;

position = [700 700 400 600];
positionNoY = [700 700 400 600];
positionNoYBar = [700 700 460 600];
positionBar = [700 700 480 600];
set(0,'DefaultFigurePaperPositionMode','auto');
colors = max(0, parula(3) - 0.1);

%==============================
% Sound speed
%==============================
Rgrid.plot_soundSpeed();
title('');
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', positionBar);
colorbar();
box on;
saveas(gcf, 'Example30_C', 'png'); 
saveas(gcf, 'Example30_C.fig'); 

%==============================
% Initial pressure
%==============================
figure;
axis([0 Rgrid.xAxis(end) 0 Rgrid.yAxis(end)]);
flipGrid = Rgrid.u';
hold on;
surf(Rgrid.xAxis, Rgrid.yAxis, flipGrid, 'EdgeColor', 'none');
view(2);
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
colorbar();
box on;
saveas(gcf, 'Example30_U', 'png'); 
saveas(gcf, 'Example30_U.fig'); 

%%  %==============================
%%  % Sensor Positions
%%  %==============================
%%  figure;
%%  axis([0 Rgrid.xAxis(end) 0 Rgrid.yAxis(end)]);
%%  hold on;
%%  plot(x1(1), x1(2), 'Marker', 'o', 'Color', colors(2, :), 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 5);
%%  plot(x2(1), x2(2), 'Marker', 'o', 'Color', colors(1, :), 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 5);
%%  plot(x3(1), x3(2), 'Marker', 'o', 'Color', colors(3, :), 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 5);
%%  % Rays
%%  for j = 1:20:nRays
%%      plot(source(1).x(j, :, 1), source(1).x(j, :, 2), 'Color', colors(2, :));
%%      plot(source(2).x(j, :, 1), source(2).x(j, :, 2), 'Color', colors(1, :));
%%      plot(source(3).x(j, :, 1), source(3).x(j, :, 2), 'Color', colors(3, :));
%%  end
%%  leg = legend('Sensor 1', 'Sensor 2', 'Sensor 3');
%%  set(leg, 'Position', [0.1 0.1 0.1 0.2])
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  set(gcf, 'pos', position);
%%  box on;
%%  %saveas(gcf, 'Example30_SensorPositions', 'png');
%%  saveas(gcf, 'Example30_SensorPositions.fig');

%==============================
% Sensor 1
%==============================
figure;
axis([0 Rgrid.xAxis(end) 0 Rgrid.yAxis(end)]);
hold on;
plot(x1(1), x1(2), 'Marker', 'o', 'Color', colors(2, :), 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 5);
% Rays
for j = 1:20:nRays
    plot(source(1).x(j, :, 1), source(1).x(j, :, 2), 'Color', colors(2, :), 'Linewidth', 2);
end
xlabel('x (m)');
ylabel('y (m)');
%set(gca, 'YTick', []);
set(gcf, 'pos', position);
box on;
saveas(gcf, 'Example30_Sensor1', 'epsc');
%saveas(gcf, 'Example30_Sensor1.fig');

%==============================
% Sensor 2
%==============================
figure;
axis([0 Rgrid.xAxis(end) 0 Rgrid.yAxis(end)]);
hold on;
plot(x2(1), x2(2), 'Marker', 'o', 'Color', colors(1, :), 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 5);
% Rays
for j = 1:20:nRays
    plot(source(2).x(j, :, 1), source(2).x(j, :, 2), 'Color', colors(1, :), 'Linewidth', 2);
end
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoY);
box on;
saveas(gcf, 'Example30_Sensor2', 'epsc');
%saveas(gcf, 'Example30_Sensor2.fig');

%==============================
% Sensor 3
%==============================
figure;
axis([0 Rgrid.xAxis(end) 0 Rgrid.yAxis(end)]);
hold on;
plot(x3(1), x3(2), 'Marker', 'o', 'Color', colors(3, :), 'MarkerSize', 10, 'LineStyle', 'none', 'LineWidth', 5);
% Rays
for j = 1:20:nRays
    plot(source(3).x(j, :, 1), source(3).x(j, :, 2), 'Color', colors(3, :), 'Linewidth', 2);
end
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoY);
box on;
saveas(gcf, 'Example30_Sensor3', 'epsc');
saveas(gcf, 'Example30_Sensor3.fig');


%%  %==============================
%%  % Sources and their respective Rays
%%  %==============================
%%  figure;
%%  hold on;
%%  % Rays
%%  for j = 1:20:nRays
%%      % Sources
%%      for n = 1:nSources
%%          plot(source(n).x(j, :, 1), source(n).x(j, :, 2), 'Color', colourMapV(n));
%%          %plot(Rgrid.source(n).xPhi(1, :, j), Rgrid.source(n).xPhi(2, :, j), 'Color', colourMapV(1));
%%      end
%%  end
%%  legend('Sensor 1', 'Sensor 2', 'Sensor 3');
%%  saveas(gcf, 'Example30_Ray', 'png');
%%  saveas(gcf, 'Example30_Ray.fig');
%%  
%%  %==============================
%%  % Isotime curve
%%  %==============================
%%  nCurves = 20;
%%  xIso = Rgrid.findIsoTime(source(2), 20, 1.5e-4, 2e-4);
%%  % Figure
%%  figure;
%%  flipGrid = Rgrid.u';
%%  axis([0 x(end) 0 y(end)]);
%%  hold on;
%%  surf(Rgrid.xAxis, Rgrid.yAxis, flipGrid, 'EdgeColor', 'none');
%%  hold on;
%%  for n = 1:nCurves
%%      plot3(xIso(n, :, 1), xIso(n, :, 2), repmat(2, [1 size(xIso, 2)]), '-m');
%%  end
%%  for j = 1:20:nRays
%%      plot3(source(2).x(j, :, 1), source(2).x(j, :, 2), repmat(2, [1 Rgrid.source(2).nPoints]), '-g');
%%  end
%%  view(2);
%%  legend('Isocurves for times tMin = 1.5e-4s, tMax = 2e-4s');
%%  saveas(gcf, 'Example30_isotime', 'epsc');
%%  saveas(gcf, 'Example30_isotime.fig');

%%  %==============================
%%  % Amplitude
%%  %==============================
%%  % Figure
%%  Rgrid.plot_amplitude(2, 100, 'real');

%%  %==============================
%%  % Non-Filtered Signal
%%  %==============================
%%  % Time Signal Amplitude
%%  figure;
%%  hold on;
%%  axis([0 3e-4 -.2 1.2]);
%%  normBeam = max(source(1).aBeam);
%%  for n = 1:nSources
%%      plot(source(n).tBeam, source(n).aBeam/normBeam, 'Color', colourMapV(n));  
%%  end
%%  legend('Sensor 1', 'Sensor 2', 'Sensor 3');
%%  xlabel('t (s)');
%%  ylabel('Amplitude');
%%  grid on;
%%  box on;
%%  saveas(gcf, 'Example30_aSignal', 'epsc');
%%  saveas(gcf, 'Example30_aSignal.fig');
%%  
%==============================
% Filtered Signal
%==============================
colours = max(0, parula(3) - 0.1);
coloursKW = parula(3);
% Normalisation factor
normRT = max(source(1).aForward);
normKWave = max(sensor_data.p(1, :));

% Figure
figure;
%%% SIGNAL
axis([0 4e-4 -1.2 1.2]);
hold on;
% Plot and compare
plot(Rgrid.tForward, source(1).aForward/normRT, 'Color', colours(2, :));
plot(Rgrid.tForward, source(2).aForward/normRT, 'Color', colours(1, :));
plot(Rgrid.tForward, source(3).aForward/normRT, 'Color', colours(3, :));
plot(kgrid.t_array, sensor_data.p(1, :)/normKWave, 'Color', coloursKW(2, :), 'LineWidth', 2);
plot(kgrid.t_array, sensor_data.p(2, :)/normKWave, 'Color', coloursKW(1, :), 'LineWidth', 2);
plot(kgrid.t_array, sensor_data.p(3, :)/normKWave, 'Color', coloursKW(3, :), 'LineWidth', 2);
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 RT', 'Sensor 2 RT', 'Sensor 3 RT', 'Sensor 1 kWave', 'Sensor 2 kWave', 'Sensor 3 kWave');
grid on;
box on;
saveas(gcf, 'Example30_aSignalConv', 'epsc');
saveas(gcf, 'Example30_aSignalConv.fig');

%%  %==============================
%%  % Error
%%  %==============================
%%  colours = max(0, parula(3) - 0.1);
%%  figure;
%%  axis([0 4e-4 -.2 .2]);
%%  hold on;
%%  error(1, :) = source(1).aForward/normRT - sensor_data.p(1, :)/normKWave;
%%  error(2, :) = source(2).aForward/normRT - sensor_data.p(2, :)/normKWave;
%%  error(3, :) = source(3).aForward/normRT - sensor_data.p(3, :)/normKWave;
%%  plot(Rgrid.tForward, error(1, :), 'Color', colours(2, :), 'LineWidth', 2);
%%  plot(Rgrid.tForward, error(2, :), 'Color', colours(1, :), 'LineWidth', 2);
%%  plot(Rgrid.tForward, error(3, :), 'Color', colours(3, :), 'LineWidth', 2);
%%  xlabel('t (s)');
%%  ylabel('Error');
%%  legend('Error 1', 'Error 2', 'Error 3');
%%  grid on;
%%  box on;
%%  saveas(gcf, 'Example30_error', 'epsc');
%%  saveas(gcf, 'Example30_error.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/;
