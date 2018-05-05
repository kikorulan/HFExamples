%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 720;           % number of grid points in the x (row) direction
Ny = 1080;           % number of grid points in the y (column) direction
dx = 0.12/Nx;        % grid point spacing in the x direction [m]
dy = 0.18/Ny;        % grid point spacing in the y direction [m]
kernelSize = 33; % Multiple of 3

Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.1;
v2 = 0.1;
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
Rgrid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];
Rgrid.setUMatrix(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays & sources
nSources = 1;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
%x1 = [Nx*dx/2; Ny*dy/100]; % Source 1
%x2 = [5*Nx*dx/6; Ny*dy/100]; % Source 2
x2 = [2*Nx*dx/3; 8*Ny*dy/9]; % Point in the shadow
x3 = [Nx*dx/100; Ny*dy/3]; % Source 3


criticalAngle = 0.7460;
% Create the new sources
nRaysQ1 = 100;
Rgrid.newSource(x3, criticalAngle-5e-3, criticalAngle+5e-3, nRaysQ1, step, tauMax);
%%  nRaysQ2 = 10;
%%  Rgrid.newSource(x3, 0, criticalAngle, nRaysQ2, step, tauMax);
%%  nRaysQ3 = 100;
%%  Rgrid.newSource(x2, pi, 3*pi/2, nRaysQ3, step, tauMax);

% Sources
coord3 = Rgrid.findCoordinates(x3);
for n = 1:nSources
    Rgrid.runFMM_point(x3);
end

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = Rgrid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
colorbar;
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_C.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = Rgrid.phase';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
colorbar;
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_FMM.fig');
