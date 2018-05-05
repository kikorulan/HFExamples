% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex41_kWave_3D;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
Nz = 128;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
dz = 1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
medium.sound_speed = 1500;
c = medium.sound_speed*ones(Nx, Ny, Nz);
% compute time
dt = 1e-7;
tMax = 1e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
radi = 10;
u0 = zeros(Nx, Ny, Nz);
%u0(floor(Nx/2), floor(Ny/2), floor(Nz/2)) = 1;
u0 = addSphere(u0,   floor(Nx/2), floor(Ny/2), floor(Nz/2), radi, 1);
u0 = addSphere(u0,   floor(Nx/2), floor(Ny/4), floor(Nz/2), radi, 2);
u0 = addSphere(u0, floor(3*Nx/4), floor(Ny/2), floor(Nz/2), radi, 3);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);
source.p0 = max(0, source.p0);
%==============================
% save data
%==============================
u0_matrix = cube2matrix(source.p0);
c_matrix = cube2matrix(c);
dlmwrite('/cs/academic/phd3/frullan/Documents/C++/HighFreq_3DRT/Build/bin/initialPressure_homogeneous.dat', u0_matrix, 'delimiter', ' ');
dlmwrite('/cs/academic/phd3/frullan/Documents/C++/HighFreq_3DRT/Build/bin/soundSpeed_homogeneous.dat',       c_matrix, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, floor(Ny/2), floor(Nz/2)) = 1;

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

save sensor_data_3balls.mat kgrid medium source sensor_data;
%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex41_kWave_3D;
load sensor_data_3balls;
%==============================
% Plot sensor measurements
%==============================
figure;
plot(kgrid.t_array, sensor_data, 'Color', 'r');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
