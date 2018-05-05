% EXAMPLE29_KWAVE computes the impulse response for 
% the initial value problem given homogeneous domain

% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex29_deltaX;

close all;
%clear all;

%=========================================================================
% SIMULATION
%=========================================================================

% create the computational grid
Nx = 240;           % number of grid points in the x (row) direction
Ny = 720;           % number of grid points in the y (column) direction
dx = 1e-4;    	% grid point spacing in the x direction [m]
dy = 1e-4;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

%===============================
% Properties of the propagation medium
%===============================
% Reference sound speed
c0 = 1500;
medium.sound_speed = c0;

%===============================
% Time Array
%===============================
% create the time array
dt = 1e-8;
tMax = 6e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium_angle1.sound_speed);

%===============================
% Mask
%===============================
% define N sensors
nSensors = 18;
sensor.mask = zeros(Nx, Ny);
for i = 1:nSensors
    sensor.mask(Nx/2, (i+1)*Ny/(nSensors+2)) = 1; % First sensor is the source
end
% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

%===============================
% Initial Pressure
%===============================
% define a single source point
source.p0 = zeros(Nx, Ny);
source.p0(Nx/2, Ny/(nSensors+2)) = 1;

%===============================
% Run the simulation
%===============================
sensorKWAVE = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);

%=========================================================================
% VISUALISATION
%=========================================================================

%=====================
% Signal at the source
%=====================
figure;
hold on;
plot(kgrid.t_array, sensorKWAVE.p(1, :), '-b');
title('Signal at the source');
legend('Source');

%=====================
% Signal at the sensor
%=====================
figure;
hold on;
for i = 1:nSensors
    plot(kgrid.t_array, sensorKWAVE.p(i, :), '-r');
end
title('Sensor');
saveas(gcf, 'Example29_kWave_sensor.fig');

save sensorKWAVE.mat sensorKWAVE kgrid sensor nSensors;

%=========================================================================
% LEAST-SQUARES MINIMISATION
%=========================================================================
%load sensorKWAVE.mat;

% Compute envelope
for i = 1:nSensors
    xEnv(i) = max(sensorKWAVE.p(i, :));
    tEnv(i) = kgrid.t_array(sensorKWAVE.p(i, :) == xEnv(i));
end

% Least squares fitting
resFun = @(x) x(1)*sqrt(x(2)./(x(2)+tEnv))-xEnv;
lb = [1e-4 1e-9];
ub = [1e3 1e-3];
x0 = [200; 1e-6];

optLS.StepTolerance = 1e-10; % Options for LS
x = lsqnonlin(resFun,x0, lb, ub, optLS);
fun = @(t) x(1)*sqrt(x(2)./(x(2)+t));

% Plot result
figure;
axis([0 tMax 0 0.05]);
hold on;
plot(tEnv, xEnv, 'ro');
plot(tEnv, fun(tEnv), '-b');
legend('kWave', 'LS fit');

% Save results
deltaX = x(2);
save deltaX.mat deltaX;
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
