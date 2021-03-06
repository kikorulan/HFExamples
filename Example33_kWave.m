% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;

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
v1 = -0.1;
v2 = 0.1;
kernelSize = 30;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
M2 = v1*c0*ones(dimX, floor(dimY/2));
M3 = v2*c0*ones(dimX, floor(dimY/2));
c = M1 + [M2 M3]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.density = 1;
    
% compute time
dt = 2.5e-7;
tMax = 2e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = zeros(Nx, Ny);
u0 = makeCircle(u0, 30, 88, 10, 1);
u0 = makeCircle(u0, 30, 168, 10, 1);
u0 = makeCircle(u0, 60, 128, 10, 1);
sourceKW.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
sourceKW.p0 = smooth(kgrid, sourceKW.p0, true);


%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
%sensor.mask(1, :) = 1;
sensor.mask(1, :) = 1;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, sourceKW, sensor, input_args{:});

save sensor_data.mat kgrid sensor sourceKW medium sensor_data input_args u0;
% reset the initial pressure
source.p0 = 0;
% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

%==============================
% Time reversal
%==============================
%%  clear p0_recon;
%%  for n = 1:256
%%      sensor.mask(:) = 0;
%%      sensor.mask(1, n) = 1;
%%      sensor.time_reversal_boundary_data = sensor_data(n, :);
%%      p0_recon{n} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  end
%%  sensor.mask(1, :) = 1;
%%  sensor.time_reversal_boundary_data = sensor_data;
%%  p0_recon{257} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  save recon_data_sensors.mat p0_recon kgrid;

%==============================
% Adjoint
%==============================
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(1, :) = 1;
source_adjoint.p = fliplr(sensor_data);
p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});

%==============================
% Subsample 128 sensors
%==============================
source_adjoint_128.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint_128.p_mask(1, 1:2:end) = 1;
source_adjoint_128.p = fliplr(sensor_data(1:2:end, :));
p0_recon_adjoint_128 = kspaceFirstOrder2D(kgrid, medium, source_adjoint_128, sensor_adjoint, input_args{:});

%==============================
% Subsample 64 sensors
%==============================
source_adjoint_64.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint_64.p_mask(1, 2:4:end) = 1;
source_adjoint_64.p = fliplr(sensor_data(2:4:end, :));
p0_recon_adjoint_64 = kspaceFirstOrder2D(kgrid, medium, source_adjoint_64, sensor_adjoint, input_args{:});

%==============================
% Subsample 32 sensors
%==============================
source_adjoint_32.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint_32.p_mask(1, 4:8:end) = 1;
source_adjoint_32.p = fliplr(sensor_data(4:8:end, :));
p0_recon_adjoint_32 = kspaceFirstOrder2D(kgrid, medium, source_adjoint_32, sensor_adjoint, input_args{:});

%==============================
% Subsample 16 sensors
%==============================
source_adjoint_16.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint_16.p_mask(1, 8:16:end) = 1;
source_adjoint_16.p = fliplr(sensor_data(8:16:end, :));
p0_recon_adjoint_16 = kspaceFirstOrder2D(kgrid, medium, source_adjoint_16, sensor_adjoint, input_args{:});

% Save data
save recon_data_adjoint.mat kgrid medium source_adjoint sensor_adjoint p0_recon_adjoint p0_recon_adjoint_128 p0_recon_adjoint_64 p0_recon_adjoint_32 p0_recon_adjoint_16 input_args;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;
load sensor_data.mat;
%load recon_data.mat;
load recon_data_adjoint.mat;

%==============================
% Initial Pressure and Distr.
%==============================
figure;
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, fliplr(source.p0)' + fliplr(sensor.mask)', [-1 1]);
colormap(getColorMap);
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
saveas(gcf, 'Example33_U_kWave.fig');

%==============================
% Sound Speed
%==============================
figure;
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
surf(X, Y, medium.sound_speed', 'EdgeColor', 'none');
view(2);
saveas(gcf, 'Example33_C.fig'); 

%==============================
% Reconstruction Time Reversal
%==============================
p0_recon_norm = max(0, p0_recon)/max(p0_recon(:));
figure;
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
surf(X, Y, p0_recon_norm', 'EdgeColor', 'none');
view(2);
title('Reconstruction - Time reversal');

%==============================
% Reconstruction Adjoint
%==============================
%p0_recon_adjoint_norm = max(0, p0_recon_adjoint.p_final)/max(p0_recon_adjoint.p_final(:));
p0_recon_adjoint_norm = max(0, p0_recon_adjoint_16.p_final)/max(p0_recon_adjoint_16.p_final(:));
figure;
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
surf(X, Y, p0_recon_adjoint_norm', 'EdgeColor', 'none');
view(2);
title('Reconstruction - Adjoint');
saveas(gcf, 'Example33_kWave_recon_adjoint.fig'); 

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
