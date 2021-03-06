% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex85_3D_veins_subsampled;
%cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled;

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
cMatrix = importdata('input_data/sound_speed.dat', ' ', 0);
c = matrix2cube(cMatrix, Nz);
% Load sound speed
medium.sound_speed = c;
medium.density = 1;

% compute time
dt = 1.667e-8;
Nt = 486;
tMax = dt*(Nt-1);
kgrid.t_array = 0:dt:tMax;

% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, u0, true);
source.p0 = max(0, source.p0);
initial_pressure_veins_smooth = source.p0;
u0_smooth = cube2matrix(initial_pressure_veins_smooth);
dlmwrite('input_data/initial_pressure_veins_80x240x240_smooth.dat', u0_smooth, 'delimiter', ' ');
save input_data/initial_pressure_veins_smooth initial_pressure_veins_smooth;
plot_projection(initial_pressure_veins_smooth, dx);

%=========================================================================
% SIMULATION
%=========================================================================
%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
numberSensorsArray = 60;
xArray = round(1:(Nx-1)/(numberSensorsArray-1):Nx);
yArray = round(1:(Ny-1)/(numberSensorsArray-1):Ny);
zArray = round(1:(Nz-1)/(numberSensorsArray-1):Nz);
% YZ faces
[Y, Z] = meshgrid(yArray, zArray);
sensorYZ1 = 1    + Nx*(Y(:)-1) + Nx*Ny*(Z(:)-1);
sensor.mask(sensorYZ1) = 1;
% Number of sensors
numberSensors = sum(sensor.mask(:))

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
save input_data/sensor_data_veins_3600sensors.mat kgrid medium source sensor input_args;
% Save to disk
filename = 'input_data/Example85_forward_input_3600sensors.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);

% Call C++ code
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
system('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/kspaceFirstOrder3D-OMP -i input_data/Example85_forward_input_3600sensors.h5 -o output_data/Example85_forward_output_3600sensors.h5');


%==============================
% Adjoint
%==============================
% Read results
sensor_data = h5read('output_data/Example85_forward_output_3600sensors.h5', '/p');
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);

% Consider all sensors
source_adjoint.p_mask = sensor.mask;
source_adjoint.p = fliplr(sensor_data);
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example85_adjoint_input_3600sensors.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example85_adjoint_input_3600sensors.h5 -o output_data/Example85_adjoint_output_3600sensors.h5 --p_final');

