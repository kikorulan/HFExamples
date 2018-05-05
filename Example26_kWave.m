% Homogeneous Propagation Medium Example

%clear all;
%load gridRT_time.mat;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
%===================================
% Grid definition
%===================================
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
load cMatrix.mat;
medium.sound_speed = cMatrix;
medium.density = 1;

% compute time
kgrid.t_array = 0:grid.dt:(size(grid.source(1).pixelAPropagation, 3)-1)*grid.dt;

%===================================
% Sensor
%===================================
% Sensor mask
sensor.mask = ones(Nx, Ny);
% run the simulation
sensor.record = {'p', 'p_final'};

%===================================
% Build the source
%===================================
source.p0 = zeros(Nx, Ny);
source.p0(121, 73) = 1;

%% Run the simulation
sensorData = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% SAVE RESULTS
% =========================================================================

% Reshape
gridKWave.signal = reshape(sensorData.p, Nx, Ny, []);
gridKWave.time = kgrid.t_array;
save gridKWave.mat gridKWave;

% Scroll View
%scrollView(permute(fliplr(gridKWave), [2 1 3]), 3, [-1e-2 1e-2]);
