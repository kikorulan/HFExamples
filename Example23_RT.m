%======================================================================
%% Example for gridRT class
%======================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex23_snell;

close all;
%clear all;

%load('workspace');
% Measure computational time
tic;
start_time = clock;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 720;           % number of grid points in the x (row) direction
Ny = 1080;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
kernelSize = 33; % Multiple of 3

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

%========================================
% Ray Shooting
%========================================

% Number of rays & sources
nRays = 100;
nSources = 3;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/c0;
tStep = min(grid.dx, grid.dy)/c0/3;

tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*c0;
tauStep = min(grid.dx, grid.dy)*c0/3;

% Sources locations
x1(1, 1, 1) = grid.Nx*grid.dx/2; 
x1(1, 1, 2) = grid.Ny*grid.dy/100; 
x2(1, 1, 1) = 5*grid.Nx*grid.dx/6; 
x2(1, 1, 2) = grid.Ny*grid.dy/100; 
x3(1, 1, 1) = grid.Nx*grid.dx/100; 
x3(1, 1, 2) = grid.Ny*grid.dy/3; 

% Create the new sources
grid.setDeltaX(grid.dx/100);
grid.newSource(x1, pi/3, 2*pi/3, nRays, tStep, tMax);
grid.newSource(x1, pi/3, 2*pi/3, nRays, tauStep, tauMax);
%grid.newSource(x2, pi/2, 5*pi/6, nRays, step, tauMax);
%grid.newSource(x3, 0, 3*pi/8, nRays, step, tauMax);

% Sources
t1 = clock;
disp(strcat('Source ', int2str(1)));
grid.computeHamil(1);
t2 = clock;
disp(['  Computation time t-Hamiltonian ' num2str(etime(t2, t1))]);

disp(strcat('Source ', int2str(2)));
%grid.computeTrajecSnell(2);
grid.computeTrajec(2);
t3 = clock;
disp(['  Computation time tau-Hamiltonian ' num2str(etime(t3, t2))]);

%==================================================================================
% Plot results
%==================================================================================
close all;
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex23_snell;

x = 0:grid.dx:(grid.Nx-1)*grid.dx;
y = 0:grid.dy:(grid.Ny-1)*grid.dy;
  
%========================================
% Sound speed
%========================================
figure;
flipGrid = grid.c';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
colorbar();
title('Sound Speed');
saveas(gcf, 'Example23_C', 'epsc'); 

%========================================
% Sources and their respective Rays
%========================================
grid.plotRays([1 2], 10);
title('Trajectories');

grid.plotRays(1, 100);
title('Trajectories t-Hamiltonian');
saveas(gcf, 'Example23_Ray', 'epsc');

grid.plotRays(2, 10);
title('Trajectories tau-Hamiltonian');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

