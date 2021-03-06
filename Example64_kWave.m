% EXAMPLE27B_KWAVE_ANGLE computes the impulse response for 
% the initial value problem given an heterogeneous domain 
% with an interface forming an angle theta with respect to the
% horizontal

% Copyright (C) 2017 Kiko Rul·lan, Marta M. Betcke

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex64_GB;
clear all;
close all;

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
kernelSize = 33;
K = ones(kernelSize);
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
X = -dimX/2+1:dimX/2;
Y = -dimY/2+1:dimY/2;
[Y, X] = meshgrid(Y, X); % generate meshgrid

% Reference sound speed
c0 = 1500;
inc = 0.1;

% Medium 2 - Small Angle
theta = 3*pi/8;
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
kgrid.t_array = 0:dt:5e-5;
%[kgrid.t_array, dt] = makeTime(kgrid, medium_angle1.sound_speed);

%===============================
% Mask
%===============================
% define 4 sensor points
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2, 2*Ny/10) = 1;
sensor.mask(Nx/2, Ny - 3*Ny/10) = 1;
sensor.mask(Nx/2, Ny - 2*Ny/10) = 1;
sensor.mask(Nx/2, Ny - Ny/10) = 1;
% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

%===============================
% Initial Pressure
%===============================
% define a single source point
source.p0 = zeros(Nx, Ny);
source.p0 = addCircle(source.p0, floor(Nx/2), floor(Ny/10), 20, 1);
source.p0 = smooth(kgrid, source.p0, true);

%===============================
% Run the simulation
%===============================
sensorKWAVE = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);

%=========================================================================
% VISUALISATION
%=========================================================================
% Change folder
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex64_GB;

%=====================
% Signal at the sensor
%=====================

% Medium 2
figure;
hold on;
plot(kgrid.t_array, sensorKWAVE.p(1, :), '-r');
plot(kgrid.t_array, sensorKWAVE.p(2, :), '-b');
plot(kgrid.t_array, sensorKWAVE.p(3, :), '-g');
plot(kgrid.t_array, sensorKWAVE.p(4, :), '-m');
legend('Sensor 1', 'Sensor 2', 'Sensor 3', 'Sensor 4');
title('Sensor - Medium');
%saveas(gcf, 'Example27b_Sensor_M2.fig');


save sensorKWAVE.mat sensorKWAVE kgrid source sensor;
