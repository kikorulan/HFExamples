% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex48_time;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
v1 = -0.3;
v2 = 0.1;
kernelSize = 30;
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
dt = 2e-7;
tMax = 2.5e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

%==============================
% Build initial pressure
%==============================
xAxis = 1:length(kgrid.t_array);
position = xAxis(end)/5;
% define a single source point
source_low.p_mask = zeros(Nx, Ny);
source_mid.p_mask = zeros(Nx, Ny);
source_high.p_mask = zeros(Nx, Ny);
source_low.p_mask(floor(Nx/2),     floor(Ny/4)) = 1;
source_mid.p_mask(floor(Nx/2), floor(31*Ny/50)) = 1;
source_high.p_mask(floor(Nx/2),  floor(4*Ny/5)) = 1;
% define a time varying gaussian source
stdDev = 40;
source.p = exp(-(xAxis-position).^2/stdDev^2);
% filter the source to remove high frequencies not supported by the gri
source.p = filterTimeSeries(kgrid, medium, source.p);
source_low.p = source.p;
source_mid.p = source.p;
source_high.p = source.p;

%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = ones(Nx, Ny);
sensor.record = {'p', 'p_final'};


% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data_low = kspaceFirstOrder2D(kgrid, medium, source_low, sensor, input_args{:});
sensor_data_mid = kspaceFirstOrder2D(kgrid, medium, source_mid, sensor, input_args{:});
sensor_data_high = kspaceFirstOrder2D(kgrid, medium, source_high, sensor, input_args{:});

save sensor_data.mat kgrid sensor source_low source_mid source_high medium c0 dt sensor_data_low sensor_data_mid sensor_data_high input_args -v7.3;
%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex48_time;
load sensor_data.mat;

%==============================
% Initial time source
%==============================
figure;
plot(kgrid.t_array, source.p, 'Color', 'r');
xlabel('t (s)');

%==============================
% Pressure
%==============================
figure;
plot(kgrid.t_array, sensor_data_low.p(64, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, sensor_data_mid.p(64, :), 'Color', 'g');
plot(kgrid.t_array, sensor_data_high.p(64, :), 'Color', 'b');
legend('TS low', 'TS mid', 'TS high');
xlabel('t (s)');
ylabel('amplitude');

%==============================
% Time propagation
%==============================
sensor_data_low_reshape  = permute(fliplr(reshape(sensor_data_low.p, Nx, Ny, [])), [2 1 3]);
sensor_data_mid_reshape  = permute(fliplr(reshape(sensor_data_mid.p, Nx, Ny, [])), [2 1 3]);
sensor_data_high_reshape = permute(fliplr(reshape(sensor_data_high.p, Nx, Ny, [])), [2 1 3]);
scrollView(sensor_data_low_reshape, 3, [-0.02 0.02]);
scrollView(sensor_data_mid_reshape, 3, [-0.02 0.02]);
scrollView(sensor_data_high_reshape, 3, [-0.02 0.02]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex48_time;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
