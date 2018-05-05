% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex55_RT_3D_half;

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
% Build domain
c0 = 1500;
v1 = -0.1;
v2 = 0.1;
kernelSize = 33;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
M2 = v2*c0*ones(dimX, floor(dimY/2));
M3 = 0*M2;
c = [M1 + [M2 M3]];
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
cMatrix = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
c = repmat(cMatrix, 1, 1, Nz);
medium.sound_speed = c;
medium.density = 1;

% compute time
dt = 1e-7;
tMax = 1.5e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
radi = 10;
u0 = zeros(Nx, Ny, Nz);
u0 = addSphere(u0,   floor(Nx/2),   floor(Ny/4), floor(Nz/2), radi, 1);
u0 = addSphere(u0,   floor(Nx/2),   floor(3*Ny/4), floor(Nz/2), radi, 1);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);
source.p0 = max(0, source.p0);
%==============================
% save data
%==============================
u0_matrix = cube2matrix(source.p0);
c_matrix = cube2matrix(c);
dlmwrite('input_data/initialPressure_2balls.dat', u0_matrix, 'delimiter', ' ');
dlmwrite('input_data/soundSpeed_half.dat',        c_matrix, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
numberSensorsArray = 3;
xArray = round(1:(Nx-1)/(numberSensorsArray-1):Nx);
yArray = round(1:(Ny-1)/(numberSensorsArray-1):Ny);
zArray = round(1:(Nz-1)/(numberSensorsArray-1):Nz);
% % XY faces
% [X, Y] = meshgrid(xArray, yArray);
% sensorXY1 = X(:) + Nx*(Y(:)-1) + Nx*Ny*   (1-1);
% sensorXY2 = X(:) + Nx*(Y(:)-1) + Nx*Ny*  (Nz-1);
% sensor.mask(sensorXY1) = 1;
% sensor.mask(sensorXY2) = 1;
% % XZ faces
[X, Z] = meshgrid(xArray, zArray);
sensorXZ1 = X(:) + Nx*   (1-1) + Nx*Ny*(Z(:)-1);
sensorXZ2 = X(:) + Nx*  (Ny-1) + Nx*Ny*(Z(:)-1);
sensor.mask(sensorXZ1) = 1;
sensor.mask(sensorXZ2) = 1;
% % YZ faces
% [Y, Z] = meshgrid(yArray, zArray);
% sensorYZ1 = 1    + Nx*(Y(:)-1) + Nx*Ny*(Z(:)-1);
% sensorYZ2 = Nx   + Nx*(Y(:)-1) + Nx*Ny*(Z(:)-1);
% sensor.mask(sensorYZ1) = 1;
% sensor.mask(sensorYZ2) = 1;
% Number of sensors
numberSensors = sum(sensor.mask(:))

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

% Save to disk
filename = 'input_data/Example55_forward_input.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);

% Call C++ code
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/GCC-5.4/lib64';
system('../kspaceFirstOrder3D-OMP -i input_data/Example55_forward_input.h5 -o output_data/Example55_forward_output.h5');
save input_data/sensor_data_2balls.mat kgrid medium source sensor input_args;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex55_RT_3D_half;

% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;

%==============================
% Sound speed
%==============================
figure;
surf(x_axis, y_axis, c(:, :, 64), 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('Sound Speed');
saveas(gcf, 'Example55_C', 'png');
saveas(gcf, 'Example55_C.fig');

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
saveas(gcf, 'Example55_U', 'png');
saveas(gcf, 'Example55_U.fig');

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('output_data/Example55_forward_output.h5', '/p');
figure;
imagesc(sensor_data);
box on;


cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
