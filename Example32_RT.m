%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex32_ampl;

close all;
%clear all;

load sensor_data.mat;
 
run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
grid.setCMatrix(c0*ones(Nx, Ny));

% Set initial pressure
grid.setUMatrix(source.p0);

%========================================
% Impulse Response
%========================================
% Set time
dt = min(grid.dx, grid.dy)/c0/2;
tMax = 1.6e-5;
grid.setTime(dt, tMax);
% Compute impulse response
grid.impulseResponse2D('IV');

save gridRT_impulse.mat grid c0 dt;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;
% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 800;% 800

% Parametrisation
tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*c0;
tauStep = min(grid.dx, grid.dy)*c0/2;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/c0;
tStep = dt;

vectorSensors = 40:8:224;
% Sources locations
clear x;
nSensors = 24;
for n = 1:nSensors
    x{n} = cat(3, (96-1)*grid.dx, (vectorSensors(n)-1)*grid.dy);
end


% Sources
%%  for n = 1:nSensors
%%      % T parametrisation
%%      grid.newSource(x{n}, pi/2, 3*pi/2, nRays, tStep, tMax);
%%      grid.computeHamil(n);
%%  end


grid.computeForwardParallel(x, pi/2, 3*pi/2, nRays, tStep, tMax, false);

% Save results
save gridRT.mat grid x -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex32_ampl;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex32_ampl;

load gridRT.mat;
%load sensor_data.mat;

%==============================
% Initial Pressure
%==============================
figure;
X = 0:grid.dx:(grid.Nx-1)*grid.dx;
Y = 0:grid.dy:(grid.Ny-1)*grid.dy;
surf(X, Y, grid.u', 'EdgeColor', 'none');
view(2);
saveas(gcf, 'Example32_U.fig'); 

%==============================
% Ray trajectories
%==============================
grid.plotRays(1, 10);
grid.plotRays(nSensors, 10);

%==============================
% Time Signals
%==============================
% Normalisation factor
normRT = max(grid.source(1).aForward);
normKWave = max(sensor_data(1, :));

% Figure
figure;
hold on;
%axis([0 1.2e-4 -0.5 0.5]);
for n = 1:nSensors
    % Plot and compare
    plot(grid.tForward, grid.source(n).aForward/normRT, 'Color', 'g');
    plot(kgrid.t_array, sensor_data(n, :)/normKWave, 'Color', 'm');
end

legend('RT signals', 'kWave signals');
saveas(gcf, 'Example32_aSignalConv', 'epsc');
saveas(gcf, 'Example32_aSignalConv.fig');



%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

