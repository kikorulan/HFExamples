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
%%  medium.density = 1;
    
% compute time
dt = 1e-8;
tMax = 2e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = zeros(Nx, Ny, Nz);
u0(floor(6*Nx/10), floor(Ny/2), floor(Nz/2)) = 1;
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);


%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(floor(4*Nx/10), floor(Ny/2), floor(Nx/2)) = 1;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation


%% Save to disk
filename = 'Example41_input.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
%%  sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex41_kWave_3D;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

%======================
% Plot filter
%======================
%%  sensor_data = h5read('Example41_output.h5', '/p');
%%  figure;
%%  plot(sensor_data);
