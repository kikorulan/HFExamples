% Homogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex39_adjoint;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex38_gaussBeams;

clear all;

%=========================================================================
% SIMULATION
%=========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);
% define the properties of the propagation medium
medium.sound_speed = 1500;
medium.density = 1;

% compute time
dt = 3e-7;
tMax = 4e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = zeros(Nx, Ny);
u0 = makeCircle(u0, 30, 128, 10, 1);
u0 = makeCircle(u0, 60, 128, 10, 1);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);

% input args
input_args = {'PMLInside', false, 'Smooth', false};

%==============================
% Forward propagation
%==============================

% Define the sensors
nSensors = 256;
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;
% run the simulation
sensor.record = {'p', 'p_final'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

save sensor_data_forward.mat sensor sensor_data kgrid source medium input_args;

%==============================
% Inverse propagation
%==============================
clear all;
load sensor_data_forward.mat;
% Time reversal
source.p0 = 0;
sensor.time_reversal_boundary_data = sensor_data.p;
p0_recon_time_reversal = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
% Dirichlet data
sensor_dirichlet.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_dirichlet.record = {'p_final'};
source_dirichlet.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_dirichlet.p_mask(1, :) = 1;
source_dirichlet.p = fliplr(sensor_data.p);
source_dirichlet.p_mode = 'dirichlet';
p0_recon_dirichlet = kspaceFirstOrder2D(kgrid, medium, source_dirichlet, sensor_dirichlet, input_args{:});

save recon_data.mat p0_recon_time_reversal p0_recon_dirichlet;
%=========================================================================
% VISUALISATION
%=========================================================================
load sensor_data_forward.mat;
load recon_data.mat;
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex39_adjoint;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex39_adjoint;

%==============================
% Sensor data
%==============================
colourMapV = cool(256);
figure;
hold on;
for n = 1:256
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n, :));
end
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, 'Example39_kWave_signal.fig');

%==============================
% Initial pressure
%==============================
figure;
surf(kgrid.x_vec, kgrid.y_vec, source.p0', 'EdgeColor', 'none');
view(2); 

%==============================
% Reconstruction time reversal
%==============================
figure;
surf(kgrid.x_vec, kgrid.y_vec, p0_recon_time_reversal', 'EdgeColor', 'none');
view(2); 

%==============================
% Reconstruction dirichlet
%==============================
figure;
surf(kgrid.x_vec, kgrid.y_vec, p0_recon_dirichlet.p_final', 'EdgeColor', 'none');
view(2); 

%==============================
% Difference TR - D
%==============================
error = p0_recon_time_reversal - p0_recon_dirichlet.p_final;
figure;
surf(kgrid.x_vec, kgrid.y_vec, error', 'EdgeColor', 'none');
view(2); 


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
