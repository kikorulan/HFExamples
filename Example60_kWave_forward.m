% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex60_3D_4balls;

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
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
dz = 1e-4;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Load sound speed
cMatrix = importdata('input_data/sound_speed.dat', ' ', 0);
c = matrix2cube(cMatrix, Nz);
medium.sound_speed = c;
medium.density = 1;

% compute time
dt = 1e-8;
tMax = 1.5e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_4balls.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, u0, true);
source.p0 = max(0, source.p0);

%=========================================================================
% SIMULATION
%=========================================================================
%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
numberSensorsArray = 5;
xArray = round(1:(Nx-1)/(numberSensorsArray-1):Nx);
yArray = round(1:(Ny-1)/(numberSensorsArray-1):Ny);
zArray = round(1:(Nz-1)/(numberSensorsArray-1):Nz);
% XY faces
[X, Y] = meshgrid(xArray, yArray);
sensorXY1 = X(:) + Nx*(Y(:)-1) + Nx*Ny*   (1-1);
sensorXY2 = X(:) + Nx*(Y(:)-1) + Nx*Ny*  (Nz-1);
sensor.mask(sensorXY1) = 1;
sensor.mask(sensorXY2) = 1;
% XZ faces
[X, Z] = meshgrid(xArray, zArray);
sensorXZ1 = X(:) + Nx*   (1-1) + Nx*Ny*(Z(:)-1);
sensorXZ2 = X(:) + Nx*  (Ny-1) + Nx*Ny*(Z(:)-1);
sensor.mask(sensorXZ1) = 1;
sensor.mask(sensorXZ2) = 1;
% XZ faces
[Y, Z] = meshgrid(yArray, zArray);
sensorYZ1 = 1    + Nx*(Y(:)-1) + Nx*Ny*(Z(:)-1);
sensorYZ2 = Nx   + Nx*(Y(:)-1) + Nx*Ny*(Z(:)-1);
sensor.mask(sensorYZ1) = 1;
sensor.mask(sensorYZ2) = 1;
% Number of sensors
numberSensors = sum(sensor.mask(:))

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Save to disk
filename = 'input_data/Example60_forward_input.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);

% Call C++ code
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/GCC-5.4/lib64';
system('../kspaceFirstOrder3D-OMP -i input_data/Example60_forward_input.h5 -o output_data/Example60_forward_output.h5');
save input_data/sensor_data_4balls.mat kgrid medium source sensor input_args;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex60_3D_4balls;

% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;

%==============================
% Sound speed
%==============================
% 2D plot
figure;
surf(x_axis, y_axis, c(:, :, 64), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('Sound Speed');
saveas(gcf, 'Example60_C', 'png');
saveas(gcf, 'Example60_C.fig');

% 3D plot
plot_cube(medium.sound_speed);
%==============================
% Initial Pressure
%==============================
figure;
surf(x_axis, y_axis, source.p0(:, :, 64), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('Initial pressure');
saveas(gcf, 'Example60_U', 'png');
saveas(gcf, 'Example60_U.fig');

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('output_data/Example60_forward_output.h5', '/p');

figure;
imagesc(sensor_data);
box on;


cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
