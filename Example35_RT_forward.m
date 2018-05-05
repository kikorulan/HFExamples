%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;

close all;
%clear all;

load sensor_data.mat;
 
run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
grid.setCMatrix(medium.sound_speed);

% Set initial pressure
grid.setUMatrix(source.p0);

%========================================
% Impulse Response
%========================================
% Set time
cMax = max(grid.c(:));
dt = min(grid.dx, grid.dy)/cMax/2;
dt = 2.5e-7;
tMax = 2e-4;
grid.setTime(dt, tMax);
% Compute impulse response
grid.impulseResponse2D('IV');

save gridRT_impulse.mat grid;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;
% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 800;% 800
nSources = 256;%256

% Parametrisation
tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*cMax;
tauStep = min(grid.dx, grid.dy)*cMax;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;
tStep = dt;

% Sources locations
clear xL xR;
for n = 1:256
    xL{n} = cat(3, 0, (n-1)*(grid.Ny-1)*grid.dy/(nSources - 1));
end

for n = 1:127
    xB{n} = cat(3, n*(grid.Nx-1)*grid.dx/(128 - 1), 0);
end

%%  % Sources
%%  for n = 1:nSources
%%      grid.newSource(x{n}, -pi/2, pi/2, nRays, tStep, tMax);
%%      grid.computeHamil(n);
%%      %grid.deleteRays(n);
%%  end

grid.computeForwardParallel(xL, -pi/2, pi/2, nRays, tStep, tMax, true);
grid.computeForwardParallel(xB, 0, pi, nRays, tStep, tMax, true);

% Save aForward data
% Assign data
aForward_RT = zeros(nSources, length(grid.source(1).aForward));
for n = 1:nSources
    aForward_RT(n, :) = grid.source(n).aForward;
end

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;

%load gridRT.mat;
%load sensor_data.mat;

%%  %==============================
%%  % Initial Pressure
%%  %==============================
%%  figure;
%%  X = 0:dx:(Nx-1)*dx;
%%  Y = 0:dy:(Ny-1)*dy;
%%  surf(X, Y, grid.u', 'EdgeColor', 'none');
%%  view(2);
%%  saveas(gcf, 'Example35_U.fig'); 
  
%%  %==============================
%%  % Sources and their respective Rays
%%  %==============================
%%  figure;
%%  hold on;
%%  axis([0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy]);
%%  for n = 1:10:nRays
%%      plot(grid.source(10).x(n, :, 1), grid.source(10).x(n, :, 2), 'Color', 'r');
%%  end
%%  title('Ray trajectories - top left Sensor');
%%  saveas(gcf, 'Example35_Rays_top.fig');
%%  
%%  figure;
%%  hold on;
%%  axis([0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy]);
%%  for n = 1:10:nRays
%%      plot(grid.source(8).x(n, :, 1), grid.source(8).x(n, :, 2), 'Color', 'b');
%%  end
%%  title('Ray trajectories - bottom left Sensor');
%%  legend('Sensor 1', 'Sensor 2', 'Sensor 3');
%%  saveas(gcf, 'Example35_Rays_bottom.fig');

%%  %==============================
%%  % Beam Signals
%%  %==============================
%%  figure;
%%  hold on;
%%  for n = 1:nSources
%%      plot(grid.source(n).tBeam, grid.source(n).aBeam);
%%  end
  
%==============================
% Time Signals
%==============================
% Norms
normRT = max(grid.source(1).aForward);
normKWave = max(sensor_data(1, :));
% Number of plots
nPlots = floor(nSources/2);
nColours = cool(nPlots);
figure;
hold on;
for n = 1:nPlots
    plot(grid.tForward, grid.source(n).aForward/normRT, 'Color', nColours(n, :));
end

for n = 1:nPlots
    plot(kgrid.t_array, sensor_data(1 + ((256-1)/(nSources-1))*(n-1), :)/normKWave, 'Color', 'r');
end

%==============================
% Time Signals
%==============================
figure;
for n = 1:nSources
    plot(grid.tForward, grid.source(n).aForward/normRT, 'Color', 'r');
    hold on;
    plot(kgrid.t_array, sensor_data(n, :)/normKWave, 'Color', 'b');
    hold off;
    %pause();
end

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

