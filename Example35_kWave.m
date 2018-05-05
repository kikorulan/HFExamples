% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;

clear all;
close all;

run colourMap;
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
kgrid.t_array = 0:dt:tMax+dt;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = zeros(Nx, Ny);
u0 = makeCircle(u0, 30, 88, 10, 1);
u0 = makeCircle(u0, 30, 168, 10, 1);
u0 = makeCircle(u0, 60, 128, 10, 1);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);

% Define the sensors
%sensor.mask = gridAux.u;
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;
sensor.mask(:, 1) = 1;
% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

save sensor_data.mat sensor_data kgrid source sensor medium u0;

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

save recon_data.mat p0_recon kgrid;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;
load sensor_data.mat;
load recon_data.mat;


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
saveas(gcf, 'Example35_U.fig'); 

%==============================
% Initial Pressure
%==============================
p0_recon_norm = max(0, p0_recon)/max(p0_recon(:));
figure;
X = 0:dx:(Nx-1)*dx;
Y = 0:dy:(Ny-1)*dy;
surf(X, Y, p0_recon_norm', 'EdgeColor', 'none');
view(2);
title('Initial Pressure');
saveas(gcf, 'Example35_U.fig'); 

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
