% Build matrix of propagation times from each pair  of sensors
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex58_kWave_US;
cd /home/kiko/Documents/HighFreqCode/Examples/Ex58_kWave_US;

%close all;
%=========================================================================
% SIMULATION
%=========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 1e-4;          % grid point spacing in the x direction [m]
dy = 1e-4;          % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
% Build domain
c0 = 1500;
factor = 0.1;
% Build Peaks
frame = 2.5;
x = -frame:(2*frame)/(Ny-1):frame;
y = -frame:(2*frame)/(Nx-1):frame;
[X, Y] = meshgrid(x, y);
p = peaks(X, Y);
medium.sound_speed = c0*(ones(Nx, Ny) + factor*p/max(p(:)));
medium.density = 1;
% Define time vector
dt = 5e-8;
tMax = 4e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
lSignal = length(kgrid.t_array);
%===============================
% Sensor
%===============================
subsampleFactor = 4;
sensor.mask = zeros(Nx, Ny);
sensor.mask(1:subsampleFactor:end, 1) = 1;     % y = 1
sensor.mask(end:-subsampleFactor:1, end) = 1;  % y = end
sensor.mask(1, 1:subsampleFactor:end) = 1;     % x = 1
sensor.mask(end, end:-subsampleFactor:1) = 1;  % x = end

nSensors = sum(sensor.mask(:));
sensor_vec = find(sensor.mask == 1);
% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

%===============================
% Time Varying Source
%===============================
width = floor(lSignal/20);
% define a time varying sinusoidal source
subsignal = 2*width:3*width;
lSub = length(subsignal);
source.p = zeros(1, lSignal);
source.p(subsignal) = exp(-(subsignal-2.5*width).^2/((width/2).^2));
% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);
% define mask
source.p_mask = zeros(Nx, Ny);
% Input parameters
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false, 'PMLSize', 20};

% Loop over sources
prop_time_matrix = zeros(nSensors, nSensors);
for i = 1:nSensors
    disp(i)
    % Single source point
    source.p_mask(:) = 0;
    source.p_mask(sensor_vec(i)) = 1;
    % Run the simulation
    sensor_data{i} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    % Extract propagation times
    prop_time_int = zeros(nSensors, 1);
    % Loop over sensors
    for j = 1:nSensors
        sensor_data_vec = sensor_data{i}.p(j, :);
        prop_time_int(j) = find(sensor_data_vec == max(sensor_data_vec));
    end
    prop_time_int = prop_time_int - min(prop_time_int);
    prop_time = dt*prop_time_int;
    prop_time_matrix(:, i) =  prop_time;
end
%=========================================================================
% VISUALISATION
%=========================================================================
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
%==============================
% Sound Speed
%==============================
figure;
surf(X, Y, medium.sound_speed', 'EdgeColor', 'none');
colorbar();
view(2);

%==============================
% Initial signal
%==============================
figure;
plot(kgrid.t_array, source.p);

%==============================
% Sensor Data
%==============================
figure;
%surf(kgrid.t_array, 1:nSensors, sensor_data.p, 'EdgeColor', 'none');
imagesc(kgrid.t_array, 1:nSensors, sensor_data.p(1:end, :));
view(2);
