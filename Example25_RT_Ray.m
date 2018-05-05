% Compute the rays and amplitude at the given domain
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 240;           % number of grid points in the x (row) direction
Ny = 240;           % number of grid points in the y (column) direction
dx = 0.1/Nx;        % grid point spacing in the x direction [m]
dy = 0.1/Ny;        % grid point spacing in the y direction [m]
kernelSize = 33;

grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.2;
v2 = 0.2;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, floor(dimY/3));
M2 = v1*c0*ones(floor(dimX/2), floor(dimY/3));
M3 = v2*c0*ones(floor(dimX/2), floor(dimY/3));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [M1...
     [M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
grid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/4),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M2 ; M3] M4];
grid.setUMatrix(U);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of rays & sources
nRays = 2000;
nSources = 1;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
x1 = [Nx*dx/2; Ny*dy/5]; % Source 1

% Create the new sources
grid.newSource(x1, 0, 2*pi-0.01, nRays, step, tauMax);

% Compute the rays and amplitudes
disp(strcat('Source ', int2str(1)));
grid.computeSourceQ(1, grid.dx, -inf);
save gridRT_Ray.mat grid nRays nSources; % Save variables

%%%%%%%%%%%%%%%%%%%%%%%%
% Time
%%%%%%%%%%%%%%%%%%%%%%%%
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

