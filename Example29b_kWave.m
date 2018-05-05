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
inc = 0.2;
kernelSize = 33;
K = ones(kernelSize);
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
X = -dimX/2+1:dimX/2;
Y = -dimY/2+1:dimY/2;
[Y, X] = meshgrid(Y, X); % generate meshgrid
theta = 0;
cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
cMatrix = c0*(1+inc*cMatrix); 
cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.density = 1;

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
saveas(gcf, 'Example29b_kWave_sensor.fig');

save sensorKWAVE.mat sensorKWAVE kgrid sensor nSensors dt tMax;

%=========================================================================
% LEAST-SQUARES MINIMISATION
%=========================================================================
%load sensorKWAVE.mat;

% Compute envelope
for i = 1:nSensors
    xEnv(i) = max(sensorKWAVE.p(i, :));
    tEnv(i) = kgrid.t_array(sensorKWAVE.p(i, :) == xEnv(i));
end

% Save results
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
