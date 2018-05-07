% Comparison between Time Varying Source vS Initial Value Source

clear all;
close all;
% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 60;           % number of grid points in the x (row) direction
Ny = 60;           % number of grid points in the y (column) direction
dx = 1e-3/Nx;      % grid point spacing in the x direction [m]
dy = 1e-3/Ny;      % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.density = 1;

% Define time vector
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

%% Sensor definition
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, Ny/2) = 1;
% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

% ===============================
% Time Varying Source
% ===============================
% define a single source point
source_TS.p_mask = zeros(Nx, Ny);
source_TS.p_mask(end - Nx/4, Ny/2) = 1;
% define a time varying sinusoidal source
source_TS.p = zeros(1, length(kgrid.t_array));
source_TS.p(floor(length(kgrid.t_array)/5)) = 1;
% filter the source to remove high frequencies not supported by the grid
source_TS.p = filterTimeSeries(kgrid, medium, source_TS.p);

% ===============================
% Initial Value Problem
% ===============================
% define a single source point
source_IV.p0 = zeros(Nx, Ny);
source_IV.p0(end - Nx/4, Ny/2) = 1;

% Run the simulation
sensorData_TS = kspaceFirstOrder2D(kgrid, medium, source_TS, sensor);
sensorData_IV = kspaceFirstOrder2D(kgrid, medium, source_IV, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulated sensor data
figure;
hold on;
plot(kgrid.t_array, sensorData_TS.p, 'Color', 'r');
plot(kgrid.t_array, sensorData_IV.p, 'Color', 'g');
xlabel('Time');
ylabel('Signal Amplitude');
grid on;
title('Sensor Pressure Signal');
legend('TS', 'IV');
saveas(gcf, 'Example15_kWave_TSvsIV.fig'); 
