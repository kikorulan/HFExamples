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
Nx = 120;
dx = 1e-3;
Ny = 180;
dy = 1e-3;
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

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];
grid.setUMatrix(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays & sources
nRays = 2;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
x0 = [Nx*dx/2; Ny*dy/100]; % Source 1
% Create the new sources
grid.newSource(x0, 3*pi/7, 5*pi/7, nRays, step, tauMax);

% Sources
grid.computeSource(1);
grid.computeAmplitudeProximal(1);
%grid.computeAmplitudeProximal(1);

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed of sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.c';
surf(x, y, flipGrid, 'EdgeColor', 'none');
colorbar;
view(2);
title('Speed of Sound');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex19_RT_interf_loss/Example19_C.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources and their respective Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
plot(grid.source(1).x(1, :, 1), grid.source(1).x(2, :, 1), 'Color', colourMapV(1)); % Main Ray 1
plot(grid.source(1).x(1, :, 2), grid.source(1).x(2, :, 2), 'Color', colourMapV(2)); % Main Ray 2
plot(grid.source(1).xPhi(1, :, 1), grid.source(1).xPhi(2, :, 1), 'Color', colourMapV(1)); % Main Ray 1
plot(grid.source(1).xPhi(1, :, 2), grid.source(1).xPhi(2, :, 2), 'Color', colourMapV(2)); % Main Ray 2
legend('Ray 1', 'Ray 2');
title('Trajectories');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex19_RT_interf_loss/Example19_Ray.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Interface losses %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = grid.interfaceLoss(1);
figure;
hold on;
plot(grid.source(1).tau, I(:, :, 1), 'r');
plot(grid.source(1).tau, I(:, :, 2), 'b');

