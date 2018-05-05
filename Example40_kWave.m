% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex40_lipschitz;

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
rng(1);
u0 = rand(Nx, Ny, 1);
sourceKW.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
sourceKW.p0 = smooth(kgrid, sourceKW.p0, true);

%=========================================================================
% SIMULATION
%=========================================================================
% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};

%==============================
% Forward
%==============================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, sourceKW, sensor, input_args{:});
save sensor_data.mat kgrid sensor sourceKW medium sensor_data input_args u0;

%==============================
% Adjoint
%==============================
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(1, :) = 1;
source_adjoint.p = fliplr(sensor_data);
p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});
save recon_data_adjoint.mat kgrid medium source_adjoint sensor_adjoint p0_recon_adjoint input_args;

%==============================
% Loop to obtain Lipschitz
%==============================
initial_pressure.p0 = p0_recon_adjoint.p_final;
nIter = 10;
for n = 1:nIter
    forward_data = kspaceFirstOrder2D(kgrid, medium, initial_pressure, sensor, input_args{:});
    source_adjoint.p = fliplr(forward_data);
    p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});
    initial_pressure.p0 = p0_recon_adjoint.p_final;
end
%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex40_lipschitz;
load sensor_data.mat;
%load recon_data.mat;
load recon_data_adjoint.mat;

position = [700 700 320 600];
positionNoY = [700 700 300 600];
positionNoYBar = [700 700 363 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');
axisGrid = [0 gridAux.xAxis(end) 0 gridAux.yAxis(end)];


%==============================
% Initial Pressure
%==============================
figure;
surf(gridAux.xAxis, gridAux.yAxis, initial_pressure.p0', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionYBar);
%title('Initial Pressure');
saveas(gcf, 'Example40_U', 'png'); 
saveas(gcf, 'Example40_U.fig'); 


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex40_lipschitz;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
