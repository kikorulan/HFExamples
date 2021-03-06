% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex27_kWave_impulse;
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
grid.setUMatrix(source.p0);
save gridRT_IV.mat grid;

%=========================================================================================
% Ray Shooting
%=========================================================================================
%=======================================
% Create the grids for the three media
%=======================================

% Medium 1 - Horizontal Line
load gridRT_IV.mat;
grid_angle1 = grid;
clear grid;
theta = 0;
cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
cMatrix = c0*(1+inc*cMatrix); 
cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
grid_angle1.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

% Medium 2 - Small Angle
load gridRT_IV.mat;
grid_angle2 = grid;
clear grid;
theta = 3*pi/8;
cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
cMatrix = c0*(1+inc*cMatrix); 
cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
grid_angle2.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));
 
% Medium 3 - Large Angle
load gridRT_IV.mat;
grid_angle3 = grid;
theta = 4*pi/9;
cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
cMatrix = c0*(1+inc*cMatrix); 
cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
grid_angle3.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));


%========================================
% Ray Shooting Parameters
%========================================
load gridRT_IV.mat;
dt = 2e-8;
tMax = 5e-5;%sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;

% Number of rays & sources
nRays = 10000;% 800
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
% Angle 1
sensorRT_angle1_s1 = grid_angle1.newSource(sensor1, 5*pi/4, 7*pi/4, nRays, dt, tMax);
grid_angle1.computeHamil(sensorRT_angle1_s1, true);
sensorRT_angle1_s2 = grid_angle1.newSource(sensor2, 5*pi/4, 7*pi/4, nRays, dt, tMax);
grid_angle1.computeHamil(sensorRT_angle1_s2, true);
sensorRT_angle1_s3 = grid_angle1.newSource(sensor3, 5*pi/4, 7*pi/4, nRays, dt, tMax);
grid_angle1.computeHamil(sensorRT_angle1_s3, true);
sensorRT_angle1_s4 = grid_angle1.newSource(sensor4, 5*pi/4, 7*pi/4, nRays, dt, tMax);
grid_angle1.computeHamil(sensorRT_angle1_s4, true);
%%  % Angle 2
%%  sensorRT_angle2_s1 = grid_angle2.newSource(sensor1, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle2.computeHamil(sensorRT_angle2_s1, true);
%%  sensorRT_angle2_s2 = grid_angle2.newSource(sensor2, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle2.computeHamil(sensorRT_angle2_s2, true);
%%  sensorRT_angle2_s3 = grid_angle2.newSource(sensor3, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle2.computeHamil(sensorRT_angle2_s3, true);
%%  sensorRT_angle2_s4 = grid_angle2.newSource(sensor4, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle2.computeHamil(sensorRT_angle2_s4, true);
%%  % Angle 3
%%  sensorRT_angle3_s1 = grid_angle3.newSource(sensor1, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle3.computeHamil(sensorRT_angle3_s1, true);
%%  sensorRT_angle3_s2 = grid_angle3.newSource(sensor2, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle3.computeHamil(sensorRT_angle3_s2, true);
%%  sensorRT_angle3_s3 = grid_angle3.newSource(sensor3, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle3.computeHamil(sensorRT_angle3_s3, true);
%%  sensorRT_angle3_s4 = grid_angle3.newSource(sensor4, 5*pi/4, 7*pi/4, nRays, dt, tMax);
%%  grid_angle3.computeHamil(sensorRT_angle3_s4, true);
% Computations

%==================================================================================
% Save grid
%==================================================================================
%save -v7.3 sensorRT.mat sensorRT_angle1_s1 sensorRT_angle1_s2 sensorRT_angle1_s3 sensorRT_angle1_s4 ...
%                  sensorRT_angle2_s1 sensorRT_angle2_s2 sensorRT_angle2_s3 sensorRT_angle2_s4 ...
%                  sensorRT_angle3_s1 sensorRT_angle3_s2 sensorRT_angle3_s3 sensorRT_angle3_s4 ...
%                  grid_angle1 grid_angle2 grid_angle3;


%%  %=======================================
%%  % Shoot the rays
%%  %=======================================
%%  % Medium 1
%%  % Compute the rays and amplitudes
%%  disp(strcat('Medium ', int2str(1)));
%%  grid_angle1.computeSource(1);
%%  grid_angle1.amplitudeProximal(1);
%%  % Medium 2
%%  % Compute the rays and amplitudes
%%  disp(strcat('Medium ', int2str(2)));
%%  grid_angle2.computeSource(1);
%%  grid_angle2.amplitudeProximal(1);
%%  
%%  % Medium 3
%%  % Compute the rays and amplitudes
%%  nRays = 3000;
%%  grid_angle3.newSource(x1, pi/3, 2*pi/3, nRays, step, tauMax); % Create the new sources
%%  disp(strcat('Medium ', int2str(3)));
%%  grid_angle3.computeSourceQ(1, 10*grid_angle3.dx, -inf);
%%  grid_angle3.amplitudeProximal(1);
%%  
%%  %=======================================
%%  % Assign Impulse Response
%%  %=======================================
%%  % Medium 1
%%  grid_angle1.setCFilter(grid.cFilter);
%%  grid_angle1.setTFilter(grid.tFilter);
%%  grid_angle1.setDelayFilter(grid.delayFilter);
%%  grid_angle1.setFilter(grid.Filter);
%%  % Medium 2
%%  grid_angle2.setCFilter(grid.cFilter);
%%  grid_angle2.setTFilter(grid.tFilter);
%%  grid_angle2.setDelayFilter(grid.delayFilter);
%%  grid_angle2.setFilter(grid.Filter);
%%  % Medium 3
%%  grid_angle3.setCFilter(grid.cFilter);
%%  grid_angle3.setTFilter(grid.tFilter);
%%  grid_angle3.setDelayFilter(grid.delayFilter);
%%  grid_angle3.setFilter(grid.Filter);
%%  
%%  save gridRT_IV.mat grid_angle1 grid_angle2 grid_angle3 Nx Ny dx dy;  
%%  
%%  %=========================================================================================
%%  % Measure sensors
%%  %=========================================================================================
%%  load gridRT_IV;
%%  
%%  % Set time signal
%%  dt = 1e-9;
%%  signal = [1 0];
%%  tArray = [0 dt];
%%  grid_angle1.setTimeSignal(1, signal, tArray);
%%  grid_angle2.setTimeSignal(1, signal, tArray);
%%  grid_angle3.setTimeSignal(1, signal, tArray);
%%  
%%  % Define the sensor
%%  sensor1 = [Nx/2; 2*Ny/10];
%%  sensor2 = [Nx/2; Ny-3*Ny/10];
%%  sensor3 = [Nx/2; Ny-2*Ny/10];
%%  sensor4 = [Nx/2; Ny-1*Ny/10];
%%  sensorsRT = [sensor1 sensor2 sensor3 sensor4];
%%  
%%  % Compute the response
%%  sensorRT_angle1 = grid_angle1.timePropagationSensor(1, sensorsRT);
%%  sensorRT_angle2 = grid_angle2.timePropagationSensor(1, sensorsRT);
%%  sensorRT_angle3 = grid_angle3.timePropagationSensor(1, sensorsRT);

%==================================================================================
% Save grid
%==================================================================================
%%  save sensorRT.mat sensorRT_angle1 sensorRT_angle2 sensorRT_angle3 ...
%%                      grid_angle1 grid_angle2 grid_angle3 sensorsRT;

%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
