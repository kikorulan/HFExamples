% EXAMPLE27A_KWAVE_GAUSS computes the kWave impulse response for the 
% given gaussian time varying source
% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

clear all;
close all;

%==========================================================================
% SIMULATION
%==========================================================================

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

%===============================
% Define the sensors
%===============================
normSensor = 3*Nx/7;
sensor.mask = zeros(Nx, Ny);
angle = pi/2;
sensor.mask(Nx/2 + floor(normSensor*cos(angle)), Ny/2 + floor(normSensor*sin(angle))) = 1;
angle = pi/2+pi/4;
sensor.mask(Nx/2 + floor(normSensor*cos(angle)), Ny/2 + floor(normSensor*sin(angle))) = 1;
angle = pi/2+pi/3;
sensor.mask(Nx/2 + floor(normSensor*cos(angle)), Ny/2 + floor(normSensor*sin(angle))) = 1;

% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

%===============================
% Initial Value Problem
%===============================
% define a single source point
source.p0 = zeros(Nx, Ny);
source.p0(Nx/2, Ny/2) = 1;

%===============================
% Run the simulation
%===============================
sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false, 'Smooth', [false false false]);
sensorData_smooth = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);
%=========================================================================
% VISUALISATION
%=========================================================================
% Change folder
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex27_kWave_impulse;

% Plot sensors
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
colormap gray;
imagesc(kgrid.y_vec, kgrid.x_vec, sensor.mask + source.p0, [-1, 1]);
saveas(gcf, 'Example27c_SensorLocations.fig');

% Plot sensor signals
figure;
subplot(3, 1, 1)
plot(kgrid.t_array, sensorData.p(1, :), '-b');
title('Sensor 1 - non-Smooth - angle = 0');
xlabel('t (s)');
grid on;

subplot(3, 1, 2)
plot(kgrid.t_array, sensorData.p(2, :), '-r')
title('Sensor 2 - non-Smooth - angle = pi/4');
xlabel('t (s)');
grid on;

subplot(3, 1, 3)
plot(kgrid.t_array, sensorData.p(3, :), '-g');
title('Sensor 3 - non-Smooth - angle = pi/3');
xlabel('t (s)');
grid on;
saveas(gcf, 'Example27c_Sensors_nonSmooth.fig');

% Plot sensor signals
figure;
subplot(3, 1, 1)
plot(kgrid.t_array, sensorData_smooth.p(1, :), '-b');
title('Sensor 1 - Smooth - angle = 0');
xlabel('t (s)');
grid on;

subplot(3, 1, 2)
plot(kgrid.t_array, sensorData_smooth.p(2, :), '-r');
title('Sensor 2 - Smooth - angle = pi/4');
xlabel('t (s)');
grid on;

subplot(3, 1, 3)
plot(kgrid.t_array, sensorData_smooth.p(3, :), '-g');
title('Sensor 3 - Smooth - angle = pi/3');
xlabel('t (s)');
grid on;
saveas(gcf, 'Example27c_Sensors_Smooth.fig');


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
