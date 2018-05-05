% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex46_caustics;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================

% create the computational grid
Nx = 256;           % number of grid points in the x (row) direction
Ny = 512;           % number of grid points in the y (column) direction
dx = 5e-4;        % grid point spacing in the x direction [m]
dy = 5e-4;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
v1 = -0.3;
v2 = 0.1;
kernelSize = 60;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
c = M1;
c = addCircle(c, kernelSize + floor(Nx/2), kernelSize + floor(Ny/2), floor(Nx/6), c0*v1);
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.density = 1;
    
% compute time
%dt = 2e-7;
%tMax = 2e-4;
%kgrid.t_array = 0:dt:tMax;
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = zeros(Nx, Ny);
u0_low  = addCircle(u0, floor(Nx/2), floor(Ny/4), floor(Nx/8), 1);
u0_mid  = addCircle(u0, floor(Nx/2), floor(3*Ny/5), floor(Nx/8), 1);
u0_high = addCircle(u0, floor(Nx/2), floor(4*Ny/5), floor(Nx/8), 1);
% smooth the initial pressure distribution and restore the magnitude
source_low.p0 = smooth(kgrid, u0_low, true);
source_mid.p0 = smooth(kgrid, u0_mid, true);
source_high.p0 = smooth(kgrid, u0_high, true);

%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
%sensor.mask(1, :) = 1;
sensor.mask(floor(Nx/2), 1) = 1;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data_low = kspaceFirstOrder2D(kgrid, medium, source_low, sensor, input_args{:});
sensor_data_mid = kspaceFirstOrder2D(kgrid, medium, source_mid, sensor, input_args{:});
sensor_data_high = kspaceFirstOrder2D(kgrid, medium, source_high, sensor, input_args{:});

save sensor_data_oversample.mat kgrid sensor source_low source_mid source_high medium c0 dt sensor_data_low sensor_data_mid sensor_data_high input_args u0;
%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex46_caustics;
load sensor_data_oversample.mat;

%==============================
% Pressure
%==============================
figure;
plot(kgrid.t_array, sensor_data_low, 'Color', 'r');
hold on;
plot(kgrid.t_array, sensor_data_mid, 'Color', 'g');
plot(kgrid.t_array, sensor_data_high, 'Color', 'b');
legend('IP low', 'IP mid', 'IP high');
xlabel('t (s)');
ylabel('amplitude');

load sensor_data.mat;
plot(kgrid.t_array, sensor_data_low, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_mid, 'Color', [0.2 0.8 0.2], 'LineWidth', 2);
plot(kgrid.t_array, sensor_data_high, 'Color', [0.2 0.2 0.8], 'LineWidth', 2);



cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex46_caustics;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
