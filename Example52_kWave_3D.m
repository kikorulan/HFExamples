% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex52_RT_3D;

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
M1 = c0*ones(floor(dimX/2), dimY);
M2 = v1*c0*ones(floor(dimX/2), floor(dimY/2));
M3 = v2*c0*ones(floor(dimX/2), floor(dimY/2));
M4 = tril(M2);
M5 = fliplr(tril(M3));
c = [[M1 + [M4 M5]];...
     [M1 + [M2 M3]]]; 
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
u0 = addSphere(u0,   floor(Nx/2),   floor(Ny/2), floor(Nz/2), radi, 1);
u0 = addSphere(u0, floor(3*Nx/4),   floor(Ny/4), floor(Nz/2), radi, 1);
u0 = addSphere(u0, floor(3*Nx/4), floor(3*Ny/4), floor(Nz/2), radi, 1);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);
source.p0 = max(0, source.p0);
%==============================
% save data
%==============================
u0_matrix = cube2matrix(source.p0);
c_matrix = cube2matrix(c);
dlmwrite('input_data/initialPressure_3balls.dat', u0_matrix, 'delimiter', ' ');
dlmwrite('input_data/soundSpeed_10p.dat',       c_matrix, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask( 1,  1,  1) = 1;
sensor.mask(Nx,  1,  1) = 1;
sensor.mask( 1, Ny,  1) = 1;
sensor.mask(Nx, Ny,  1) = 1;
sensor.mask( 1,  1, Nz) = 1;
sensor.mask(Nx,  1, Nz) = 1;
sensor.mask( 1, Ny, Nz) = 1;
sensor.mask(Nx, Ny, Nz) = 1;
%sensor.mask(1, :, :) = 1;
%sensor.mask(128, :, :) = 1;

% set the input arguments: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
%sensor_data_3balls = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%save sensor_data_3balls_inhomogeneous.mat kgrid medium source sensor_data_3balls;

% Save to disk
filename = 'input_data/Example52_forward_input.h5';
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);

% Call C++ code
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/GCC-5.4/lib64';
system('../kspaceFirstOrder3D-OMP -i input_data/Example52_forward_input.h5 -o output_data/Example52_forward_output.h5');

%==============================
% Adjoint
%==============================
%% Read results
sensor_data_3balls = h5read('output_data/Example52_forward_output.h5', '/p');
save input_data/sensor_data.mat sensor_data_3balls;
source_adjoint.p_mask = sensor.mask;
source_adjoint.p = fliplr(sensor_data_3balls);
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example52_adjoint_input.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example52_adjoint_input.h5 -o output_data/Example52_adjoint_output.h5 --p_final');

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex52_RT_3D;

% Axis
x_axis = 0:dx:(Nx-1)*dx;
y_axis = 0:dy:(Ny-1)*dy;

%==============================
% Reconstruction
%==============================
p0_recon_PML = h5read('output_data/Example52_adjoint_output.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));

% Normalisation - kWave
inputKWave = p0_recon(:, :, 64);
normKWave = max(inputKWave(:));
inputKWave_norm = inputKWave/normKWave;

% Plot figure
figure;
surf(x_axis, y_axis, inputKWave_norm, 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('kWave recon - z = 0.063');
saveas(gcf, 'Example52_kWave_recon', 'png');
saveas(gcf, 'Example52_kWave_recon.fig');

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
saveas(gcf, 'Example52_C', 'png');
saveas(gcf, 'Example52_C.fig');

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
saveas(gcf, 'Example52_U', 'png');
saveas(gcf, 'Example52_U.fig');


cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex52_RT_3D;
