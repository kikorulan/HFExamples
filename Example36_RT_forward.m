%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;

close all;
%clear all;

load sensor_data.mat;
 
run colourMap;
%%  %========================================
%%  % Grid definition
%%  %========================================
%%  Nx = 128;           % number of grid points in the x (row) direction
%%  Ny = 256;           % number of grid points in the y (column) direction
%%  dx = 1e-3;        % grid point spacing in the x direction [m]
%%  dy = 1e-3;        % grid point spacing in the y direction [m]
%%  grid = gridRT(Nx, dx, Ny, dy);
%%  
%%  % Build domain
%%  grid.setCMatrix(medium.sound_speed);
%%  
%%  % Set initial pressure
%%  grid.setUMatrix(source.p0);
%%  
%%  %========================================
%%  % Impulse Response
%%  %========================================
%%  % Set time
    cMax = max(grid.c(:));
%%  %dt = min(grid.dx, grid.dy)/cMax/2;
    dt = 2.5e-7;
%%  tMax = 2e-4;
%%  grid.setTime(dt, tMax);
%%  % Compute impulse response
%%  grid.impulse_additive('IV');
%%  
%%  save gridRT_impulse.mat grid cMax dt;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;
% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 800;% 800
nSources = 764;%256

% Parametrisation
tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*cMax;
tauStep = min(grid.dx, grid.dy)*cMax;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;
tStep = dt;

% Sources locations
clear xL xR xB xT;
for n = 1:256
    xL{n} = cat(3, 0, (n-1)*(grid.Ny-1)*grid.dy/(256 - 1));
    xR{n} = cat(3, (grid.Nx-1)*grid.dx, (n-1)*(grid.Ny-1)*grid.dy/(256 - 1));
end

for n = 1:126
    xB{n} = cat(3, n*(grid.Nx-1)*grid.dx/(128 - 1), 0);
    xT{n} = cat(3, n*(grid.Nx-1)*grid.dx/(128 - 1), (grid.Ny-1)*grid.dy);
end

%%  % Sources
%%  for n = 1:nSources
%%      grid.newSource(x{n}, -pi/2, pi/2, nRays, tStep, tMax);
%%      grid.computeHamil(n);
%%      %grid.deleteRays(n);
%%  end

source = grid.computeForwardParallel(xL, -pi/2, pi/2, nRays, tStep, tMax, false);
source = [source grid.computeForwardParallel(xR, pi/2, 3*pi/2, nRays, tStep, tMax, false)];
source = [source grid.computeForwardParallel(xB, 0, pi, nRays, tStep, tMax, false)];
source = [source grid.computeForwardParallel(xT, pi, 2*pi, nRays, tStep, tMax, false)];

% Save aForward data
% Assign data
aForward_RT = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    aForward_RT(n, :) = source(n).aForward;
end

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;

%load gridRT.mat;
%load sensor_data.mat;

%%  %==============================
%%  % Time Signals
%%  %==============================
%%  % Norms
%%  normRT = max(grid.source(1).aForward);
%%  normKWave = max(sensor_data(1, :));
%%  % Number of plots
%%  nPlots = floor(nSources/2);
%%  nColours = cool(nPlots);
%%  figure;
%%  hold on;
%%  for n = 1:nPlots
%%      plot(grid.tForward, grid.source(n).aForward/normRT, 'Color', nColours(n, :));
%%  end
%%  
%%  for n = 1:nPlots
%%      plot(kgrid.t_array, sensor_data(1 + ((256-1)/(nSources-1))*(n-1), :)/normKWave, 'Color', 'r');
%%  end
%%  
%%  %==============================
%%  % Time Signals
%%  %==============================
%%  figure;
%%  for n = 1:nSources
%%      plot(grid.tForward, grid.source(n).aForward/normRT, 'Color', 'r');
%%      hold on;
%%      plot(kgrid.t_array, sensor_data(n, :)/normKWave, 'Color', 'b');
%%      hold off;
%%      %pause();
%%  end

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

