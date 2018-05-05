% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex26_timeStep_imp;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex26_timeStep_imp;

close all;
%clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
%=======================================
% Grid definition
%=======================================
Nx = 240;           
Ny = 360;           
dx = 1e-4;        
dy = 1e-4;        
kernelSize = 33;

grid = gridRT(Nx, dx, Ny, dy);

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
cMatrix = conv2(c, K, 'same')/kernelSize/kernelSize;
cMatrix = cMatrix(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
grid.setCMatrix(cMatrix);
save cMatrix.mat cMatrix;

%========================================
% Impulse Response
%========================================
% Set time
cMax = max(grid.c(:));
dt = min(grid.dx, grid.dy)/cMax/3;
tMax = 2e-5;
grid.setTime(dt, tMax);
% Compute impulse response
grid.impulseResponse2D('IV');

save gridRT_impulse.mat grid;
%=======================================
% Ray Shooting
%=======================================
% Number of rays & sources
nRays = 2000; % 3000
nSources = 1;
% Maximum t and step size
tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;
tStep = dt;

% Sources locations
x1 = cat(3, grid.Nx*grid.dx/2, grid.Ny*grid.dy/5); % Source 1

% Create the new sources
grid.newSource(x1, 0, 2*pi-0.01, nRays, tStep, tMax);

% Compute the rays and amplitudes
disp(strcat('Source ', int2str(1)));
grid.computeHamil(1);
grid.deleteRays(1);

save gridRT_Ray.mat grid nRays nSources -v7.3; % Save variables

%=======================================
% Sound Speed
%=======================================

% Source Matrix
coordS = grid.findCoordinates(x1);
sourceMatrix = zeros(Nx, Ny);
sourceMatrix(coordS(1)-2:coordS(1)+2, coordS(2)-2:coordS(2)+2) = inf;
sourceMatrix(coordS(1)-1:coordS(1)+1, coordS(2)-1:coordS(2)+1) = -300;

% Medium 1
figure;
sourceCMatrix = sourceMatrix + grid.c;
surf(grid.xAxis, grid.yAxis, sourceCMatrix', 'EdgeColor', 'none');
view(2);
title('Sound Speed & Source');
xlabel('x (m)');
ylabel('y (m)');
colorbar();
saveas(gcf, 'Example26_C_SensorLocations.fig');

%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/;
