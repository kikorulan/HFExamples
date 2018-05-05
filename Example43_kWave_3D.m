% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex43_RT_3D;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================

% create the computational grid
Nx = 64;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
Nz = 128;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
dz = 1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
v1 = -0.1;
v2 = 0.1;
kernelSize = 33;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, floor(dimY/2));
M2 = v1*c0*ones(dimX, floor(dimY/2));
M3 = v2*c0*ones(dimX, floor(dimY/2));
c = [[M1 + M2] [M1 + M3]]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
cMatrix = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
c = repmat(cMatrix, 1, 1, Nz);
medium.sound_speed = c;
medium.density = 1;

% compute time
dt = 1e-7;
tMax = 1e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
radi = 10;
u0 = zeros(Nx, Ny, Nz);
%u0(floor(Nx/2), floor(Ny/2), floor(Nz/2)) = 1;
u0 = addSphere(u0, floor(Nx/2),   floor(Ny/4), floor(Nz/2), radi, 1);
u0 = addSphere(u0, floor(Nx/2), floor(3*Ny/4), floor(Nz/2), radi, 1);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);
source.p0 = max(0, source.p0);
%==============================
% save data
%==============================
u0_matrix = cube2matrix(source.p0);
c_matrix = cube2matrix(c);
dlmwrite('/cs/academic/phd3/frullan/Documents/C++/HighFreq_3DRT/Build/bin/initialPressure_2balls_rectangle.dat', u0_matrix, 'delimiter', ' ');
dlmwrite('/cs/academic/phd3/frullan/Documents/C++/HighFreq_3DRT/Build/bin/soundSpeed_rectangle.dat', c_matrix, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(floor(Nx/2), 1, floor(Nz/2)) = 1;

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
%sensor_data_rectangle = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

%% Save to disk
filename = 'example_input.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);



save sensor_data_rectangle.mat kgrid medium source sensor_data_rectangle;
%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex43_RT_3D;
load sensor_data_rectangle;
%==============================
% Plot sensor measurements
%==============================
figure;
plot(kgrid.t_array, sensor_data_rectangle, 'Color', 'r');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
