% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex29_deltaX;

close all;
clear all;

load sensorKWAVE.mat;

% Measure computational time
tic;
start_time = clock;
  
%=======================================
% Grid definition
%=======================================
Nx = 240;           
Ny = 720;           
dx = 1e-4;        
dy = 1e-4;        

grid = gridRT(Nx, dx, Ny, dy);

% ===============================
% Properties of the propagation medium
% ===============================
% Reference sound speed
c0 = 1500;
grid.setCMatrix(c0*ones(Nx, Ny));

%=======================================
% Impulse response
%=======================================
% Compute impulse response
grid.impulseResponse2D(1e-8, 'IV', true);
save gridRT_IV.mat grid;

%=========================================================================================
% Ray Shooting
%=========================================================================================
load gridRT_IV.mat;
% Number of rays & sources
nRays = 1;
nSources = 1;
% Maximum tau and step size
tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;
% Sources locations
x1 = [(Nx/2-1)*dx; (Ny/(nSensors+2)-1)*dy]; % Source 1

%=======================================
% Shoot the rays                        
%=======================================
grid.setDeltaX(grid.dx/100);
grid.newSource(x1, pi/2, pi/2, nRays, step, tauMax); % Create the new sources
grid.computeSource(1);
grid.amplitudeProximal(1);

save gridRT_IV.mat grid;  

%=========================================================================================
% Measure sensors
%=========================================================================================
load gridRT_IV;
% Set time signal
dt = 1e-8;
signal = [1 0];
tArray = [0 dt];
grid.setTimeSignal(1, signal, tArray);

% Define the sensor
for i = 1:nSensors
    sensorRT(1, i) = grid.Nx/2;
    sensorRT(2, i) = (i+1)*grid.Ny/(nSensors+2);
end

% Compute the response
sensorRT_signal = grid.timePropagationSensor(1, sensorRT);

%==================================================================================
% Save grid
%==================================================================================
save sensorRT.mat sensorRT_signal sensorRT;

%=========================================================================================
% Plot results
%=========================================================================================
% Normalisation factors
normKWAVE = max(sensorKWAVE.p(1, :));
normRT = max(sensorRT_signal(1, :));

figure;
hold on;
for i = 1:nSensors
    plot(grid.timeSignal, sensorRT_signal(i, :)/normRT, '-r', 'LineWidth', 2);
    plot(kgrid.t_array, sensorKWAVE.p(i, :)/normKWAVE, '-b');
end
legend('RT signals', 'kWave signals');
saveas(gcf, 'Example29_sensorComparison_homogen.fig');

%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
