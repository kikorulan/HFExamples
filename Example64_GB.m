% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex64_GB;
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

%%  load gridRT_IV;

%=======================================
% Grid definition
%=======================================
Nx = 240;           
Ny = 720;           
dx = 1e-4;        
dy = 1e-4;        

grid = gridRT(Nx, dx, Ny, dy);

%===============================
% Properties of the propagation medium
%===============================
kernelSize = 33;
K = ones(kernelSize);
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
X = -dimX/2+1:dimX/2;
Y = -dimY/2+1:dimY/2;
[Y, X] = meshgrid(Y, X); % generate meshgrid
% Reference sound speed
c0 = 1500;
inc = 0.1;
grid.setCMatrix(c0*ones(Nx, Ny));

%=======================================
% Impulse response
%=======================================
% Compute impulse response
% Set time
cMax = max(grid.c(:));
dt = 2e-8;
tMax = 5e-5;%sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;
grid.setTime(dt, tMax);
% Compute impulse response
grid.impulse_additive('IV');

%===============================
% Initial Pressure
%===============================
load sensorKWAVE;
u0 = source.p0;
grid.setUMatrix(u0);
save gridRT_IV.mat grid;

%=========================================================================================
% Ray Shooting
%=========================================================================================
%=======================================
% Create the grids for the medium
%=======================================
load gridRT_IV.mat;
theta = 3*pi/8;
cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
cMatrix = c0*(1+inc*cMatrix); 
cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
grid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

%========================================
% Ray Shooting Parameters
%========================================
dt = 2e-8;
tMax = 5e-5;%sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;

% Number of rays & sources
nRays = 2000;% 800
%========================================
% Forward Problem
%========================================
% Measure computational time
tic;
start_time = clock;

% Sources locations
clear x;
sensor1 = cat(3, (grid.Nx-1)*grid.dx/2, 2*(grid.Ny-1)*grid.dy/10);
sensor2 = cat(3, (grid.Nx-1)*grid.dx/2, (grid.Ny-1)*grid.dy-3*(grid.Ny-1)*grid.dy/10);
sensor3 = cat(3, (grid.Nx-1)*grid.dx/2, (grid.Ny-1)*grid.dy-2*(grid.Ny-1)*grid.dy/10);
sensor4 = cat(3, (grid.Nx-1)*grid.dx/2, (grid.Ny-1)*grid.dy-1*(grid.Ny-1)*grid.dy/10);
% Ray Tracing
sensorRT_s1 = grid.newSource(sensor1, 5*pi/4, 7*pi/4, nRays, dt, tMax);
sensorRT_s2 = grid.newSource(sensor2, 5*pi/4, 7*pi/4, nRays, dt, tMax);
sensorRT_s3 = grid.newSource(sensor3, 5*pi/4, 7*pi/4, nRays, dt, tMax);
sensorRT_s4 = grid.newSource(sensor4, 5*pi/4, 7*pi/4, nRays, dt, tMax);
grid.computeHamil(sensorRT_s1, 'p');
grid.computeHamil(sensorRT_s2, 'p');
grid.computeHamil(sensorRT_s3, 'p');
grid.computeHamil(sensorRT_s4, 'p');
% Gaussian Beam
sensorGB_s1 = grid.newSource(sensor1, 5*pi/4, 7*pi/4, nRays, dt, tMax);
sensorGB_s2 = grid.newSource(sensor2, 5*pi/4, 7*pi/4, nRays, dt, tMax);
sensorGB_s3 = grid.newSource(sensor3, 5*pi/4, 7*pi/4, nRays, dt, tMax);
sensorGB_s4 = grid.newSource(sensor4, 5*pi/4, 7*pi/4, nRays, dt, tMax);
grid.computeHamil(sensorGB_s1, 'g');
grid.computeHamil(sensorGB_s2, 'g');
grid.computeHamil(sensorGB_s3, 'g');
grid.computeHamil(sensorGB_s4, 'g');

%==================================================================================
% Extract and save amplitude
%==================================================================================
amplitude{1} = sensorRT_s1.amplitude;
amplitude{2} = sensorRT_s2.amplitude;
amplitude{3} = sensorRT_s3.amplitude;
amplitude{4} = sensorRT_s4.amplitude;

amplitudeGB{1} = sensorGB_s1.amplitude;
amplitudeGB{2} = sensorGB_s2.amplitude;
amplitudeGB{3} = sensorGB_s3.amplitude;
amplitudeGB{4} = sensorGB_s4.amplitude;
save -v7.3 sensorRT.mat amplitude amplitudeGB grid;


%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
