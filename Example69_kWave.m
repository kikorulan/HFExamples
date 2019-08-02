% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex69_3D_RT_norm;

clear all;
close all;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
Nx = 80;           % number of grid points in the x (row) direction
Ny = 240;           % number of grid points in the y (column) direction
Nz = 240;           % number of grid points in the y (column) direction
dx = 5.3e-5;        % grid point spacing in the x direction [m]
dy = 5.3e-5;        % grid point spacing in the y direction [m]
dz = 5.3e-5;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% compute time
dt = 1.667e-8;
tMax = 8.0836e-06;
kgrid.t_array = 0:dt:tMax;

% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, u0, true);
source.p0 = max(0, source.p0);
initial_pressure_veins_smooth = source.p0;
initial_pressure_veins_smooth_mat = cube2matrix(initial_pressure_veins_smooth);
dlmwrite('input_data/initial_pressure_veins_80x240x240_smooth.dat', initial_pressure_veins_smooth, 'delimiter', ' ');
save input_data/initial_pressure_veins_smooth initial_pressure_veins_smooth;

%=========================================================================
% BUILD SOUND SPEED
%=========================================================================
% Set sound speed
c0 = 1580;
c = c0*ones(Nx, Ny, Nz);
% Save data
c0_matrix = cube2matrix(c);
dlmwrite('input_data/sound_speed.dat', c0_matrix, 'delimiter', ' ');
% Assign
medium.sound_speed = c;
medium.density = 1;
%=========================================================================
% SIMULATION
%=========================================================================
%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1, 1, 1) = 1;
sensor.mask(1, floor(Ny/4), floor(Nz/4)) = 1;
sensor.mask(1, floor(Ny/2), floor(Nz/2)) = 1;

% Number of sensors
numberSensors = sum(sensor.mask(:))

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
save input_data/sensor_data_veins_3sensors.mat kgrid medium source sensor input_args;
% Save to disk
filename = 'input_data/Example69_forward_input_3sensors.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);

% Call C++ code
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
system('../kspaceFirstOrder3D-OMP -i input_data/Example69_forward_input_3sensors.h5 -o output_data/Example69_forward_output_3sensors.h5');
