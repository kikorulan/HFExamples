% EXAMPLE27A_KWAVE_GAUSS computes the kWave impulse response for the 
% given gaussian time varying source
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

clear all;
close all;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 60;           % number of grid points in the x (row) direction
Ny = 60;           % number of grid points in the y (column) direction
dx = 1e-3;    	% grid point spacing in the x direction [m]
dy = 1e-3;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create the time array
dt = 3e-8;
kgrid.t_array = 0:dt:9e-5;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, Ny/2) = 1; % First sensor is the one away from the source
sensor.mask(end - Nx/4, Ny/2) = 1;
% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

% ===============================
% Time Varying Source
% ===============================
xAxis = 1:length(kgrid.t_array);
position = xAxis(end)/5;
% define a single source point
source_gauss.p_mask = zeros(Nx, Ny);
source_gauss.p_mask(end - Nx/4, Ny/2) = 1;
source_imp.p_mask = zeros(Nx, Ny);
source_imp.p_mask(end - Nx/4, Ny/2) = 1;
% define a time varying gaussian source
stdDev = 40;
source_gauss.p = exp(-(xAxis-position).^2/stdDev^2);
% define a single 1 time varying source
source_imp.p = zeros(1, length(kgrid.t_array));
source_imp.p(floor(position)) = 1;
% filter the source to remove high frequencies not supported by the gri
source_gauss.p = filterTimeSeries(kgrid, medium, source_gauss.p);
source_imp.p = filterTimeSeries(kgrid, medium, source_imp.p);

% ===============================
% Run the simulation
% ===============================
sensorData_gauss = kspaceFirstOrder2D(kgrid, medium, source_gauss, sensor);
sensorData_imp = kspaceFirstOrder2D(kgrid, medium, source_imp, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================
% Change folder
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex27_kWave_impulse;

% Plot the source
figure;
plot(kgrid.t_array, source_gauss.p, '-b');
hold on;
plot(kgrid.t_array, source_imp.p, '-r');
legend('Source - Gaussian', 'Source - Impulse');
grid on;
saveas(gcf, 'Example27a_Source.fig');

% Plot the signal at the sensor
axisY = [0 1e-4 -0.15 0.25];
figure;
subplot(2, 1, 1);
axis(axisY);
plot(kgrid.t_array, sensorData_gauss.p(1, :), '-b');
title('Signal at the Sensor - Gaussian source');
xlabel('t (s)');
grid on;

subplot(2, 1, 2);
axis(axisY);
plot(kgrid.t_array, sensorData_imp.p(1, :), '-r');
title('Signal at the Sensor - Impulse source');
xlabel('t (s)');
grid on;
saveas(gcf, 'Example27a_Sensors.fig');

% Comparison Source - Sensor
axisY = [0 1e-4 -0.15 0.25];
figure;
subplot(2, 1, 1);
axis(axisY);
plot(kgrid.t_array, sensorData_gauss.p(2, :), '-b');
title('Signal at the Source - Gaussian source');
xlabel('t (s)');
grid on;

subplot(2, 1, 2);
axis(axisY);
plot(kgrid.t_array, sensorData_gauss.p(1, :), '-b');
title('Signal at the Sensor - Gaussian source');
xlabel('t (s)');
grid on;
saveas(gcf, 'Example27a_SourceVsSensor.fig');

% Plot sensors
X = 0:dx:(Nx - 1)*dx;
Y = 0:dy:(Ny - 1)*dy;
figure;
colormap gray;
imagesc(X, Y, sensor.mask, [-1 1]);
saveas(gcf, 'Example27a_SensorLocation.fig');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;

%%  % ===============================
%%  % Obtain the filters
%%  % ===============================
%%  % Time Varying source
%%  maxF = max(sensorData.p(1, :));
%%  delay = find(sensorData.p(1, :) == maxF) - 1;
%%  flipFilter = fliplr(sensorData.p(1, 1:delay));
%%  indexIni = find(flipFilter < 1e-3*maxF, 1);
%%  flipFilter = fliplr(sensorData.p(1, delay+1:end));
%%  indexFin = find(abs(flipFilter) < 1e-2*maxF, 1);
%%  
%%  Filter = sensorData.p(1, delay-indexIni+1:end-indexFin);
%%  tArray = kgrid.t_array(1:end-indexFin-delay+indexIni);
