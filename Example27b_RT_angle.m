% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex27_kWave_impulse;
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

%%  load gridRT_IV;
%%  
%%  %=======================================
%%  % Grid definition
%%  %=======================================
%%  Nx = 240;           
%%  Ny = 720;           
%%  dx = 1e-4;        
%%  dy = 1e-4;        
%%  
%%  grid = gridRT(Nx, dx, Ny, dy);
%%  
%%  % ===============================
%%  % Properties of the propagation medium
%%  % ===============================
%%  kernelSize = 33;
%%  K = ones(kernelSize);
%%  dimX = Nx + 2*kernelSize;
%%  dimY = Ny + 2*kernelSize;
%%  X = -dimX/2+1:dimX/2;
%%  Y = -dimY/2+1:dimY/2;
%%  [Y, X] = meshgrid(Y, X); % generate meshgrid
%%  % Reference sound speed
%%  c0 = 1500;
%%  inc = 0.2;
%%  grid.setCMatrix(c0*ones(Nx, Ny));
%%  
%%  %=======================================
%%  % Impulse response
%%  %=======================================
%%  % Compute impulse response
%%   grid.impulseResponse2D(1e-9, 'IV', true);
%%  
%%  %=========================================================================================
%%  % Ray Shooting
%%  %=========================================================================================
%%  
%%  % Number of rays & sources
%%  nRays = 1000;
%%  nSources = 1;
%%  % Maximum tau and step size
%%  tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
%%  step = min(dx, dy)*c0/3;
%%  % Sources locations
%%  x1 = [(Nx/2-1)*dx; (Ny/10-1)*dy]; % Source 1
%%  
%%  
%%  save gridRT_IV.mat grid;
%%  clear grid;
%%  %=======================================
%%  % Create the grids for the three media
%%  %=======================================
%%  
%%  % Medium 1 - Horizontal Line
%%  load gridRT_IV.mat;
%%  grid_angle1 = grid;
%%  clear grid;
%%  theta = 0;
%%  cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
%%  cMatrix = c0*(1+inc*cMatrix); 
%%  cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
%%  grid_angle1.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));
%%  
%%  % Medium 2 - Small Angle
%%  load gridRT_IV.mat;
%%  grid_angle2 = grid;
%%  clear grid;
%%  theta = 3*pi/8;
%%  cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
%%  cMatrix = c0*(1+inc*cMatrix); 
%%  cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
%%  grid_angle2.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));
%%   
%%  % Medium 3 - Large Angle
%%  grid_angle3 = gridRT(Nx, dx, Ny, dy);
%%  theta = 4*pi/9;
%%  cMatrix = double((tan(theta) > Y./X & X < 0) | (tan(theta) < Y./X & X >= 0));
%%  cMatrix = c0*(1+inc*cMatrix); 
%%  cConv = conv2(cMatrix, K, 'same')/kernelSize/kernelSize;
%%  grid_angle3.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));
%%  grid_angle3.impulseResponse2D(1e-8, 'IV', true); % Compute impulse response;
%%    
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

%%  save gridRT_IV.mat grid_angle1 grid_angle2 grid_angle3 Nx Ny dx dy;  

%=========================================================================================
% Measure sensors
%=========================================================================================
load gridRT_IV;

% Set time signal
dt = 1e-9;
signal = [1 0];
tArray = [0 dt];
grid_angle1.setTimeSignal(1, signal, tArray);
grid_angle2.setTimeSignal(1, signal, tArray);
grid_angle3.setTimeSignal(1, signal, tArray);

% Define the sensor
sensor1 = [Nx/2; 2*Ny/10];
sensor2 = [Nx/2; Ny-3*Ny/10];
sensor3 = [Nx/2; Ny-2*Ny/10];
sensor4 = [Nx/2; Ny-1*Ny/10];
sensorsRT = [sensor1 sensor2 sensor3 sensor4];

% Compute the response
sensorRT_angle1 = grid_angle1.timePropagationSensor(1, sensorsRT);
sensorRT_angle2 = grid_angle2.timePropagationSensor(1, sensorsRT);
sensorRT_angle3 = grid_angle3.timePropagationSensor(1, sensorsRT);

%==================================================================================
% Save grid
%==================================================================================
save sensorRT.mat sensorRT_angle1 sensorRT_angle2 sensorRT_angle3 ...
                    grid_angle1 grid_angle2 grid_angle3 sensorsRT;

%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
