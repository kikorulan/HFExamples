%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex34_reconstruction;

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
c0 = 1500;
grid.setCMatrix(c0*ones(Nx, Ny));

% Set initial pressure
grid.setUMatrix(source.p0);


%========================================
% Impulse Response
%========================================
% Set time
dt = 3e-7;
%dt = min(grid.dx, grid.dy)/c0/2;
tMax = 2e-4;
grid.setTime(dt, tMax);
% Compute impulse response
grid.impulseResponse2D('IV');

save gridRT_impulse.mat grid c0;

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
tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*c0;
tauStep = min(grid.dx, grid.dy)*c0;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/c0;
tStep = dt;

% Sources locations
clear xL xR;
for n = 1:nSources
    xL{n} = cat(3, 0, (n-1)*(grid.Ny-1)*grid.dy/(nSources - 1));
    xR{n} = cat(3, (grid.Nx-1)*grid.dx, (n-1)*(grid.Ny-1)*grid.dy/(nSources - 1));
end

%%  % Sources
%%  for n = 1:nSources
%%      grid.newSource(x{n}, -pi/2, pi/2, nRays, tStep, tMax);
%%      grid.computeHamil(n);
%%      %grid.deleteRays(n);
%%  end

grid.computeForwardParallel(xL, -pi/2, pi/2, nRays, tStep, tMax, true);
grid.computeForwardParallel(xR, pi/2, 3*pi/2, nRays, tStep, tMax, true);

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
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex34_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex34_reconstruction;

%load gridRT.mat;
%load sensor_data.mat;

X = 0:grid.dx:(grid.Nx-1)*grid.dx;
Y = 0:grid.dy:(grid.Ny-1)*grid.dy;

%%  %==============================
%%  % Initial Pressure
%%  %==============================
%%  figure;
%%  surf(X, Y, grid.u', 'EdgeColor', 'none');
%%  view(2);
%%  saveas(gcf, 'Example34_U.fig'); 
%%  
%%  %==============================
%%  % Sources and their respective Rays
%%  %==============================
%%  for n = 1:nSources
%%      grid.plotRays(n, 10);
%%  end
%%  
%%  %==============================
%%  % Beam Signals
%%  %==============================
%%  figure;
%%  hold on;
%%  for n = 1:nSources
%%      plot(grid.source(n).tBeam, grid.source(n).aBeam);
%%  end
%%
  
%==============================
% Time Signals
%==============================
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

%%  nColours = cool(nPlots);
%%  figure;
%%  for n = 1:nPlots
%%      plot(grid.tForward, grid.source(n).aForward/normRT, 'r');
%%      hold on;
%%      plot(kgrid.t_array, sensor_data(1 + ((256-1)/(nSources-1))*(n-1), :)/normKWave, 'Color', 'b');
%%      errorRT = grid.source(n).aForward/normRT - sensor_data(1 + ((256-1)/(nSources-1))*(n-1), :)/normKWave;
%%      plot(grid.tForward, errorRT, 'g');
%%      hold off;
%%      legend('RT', 'kWave', 'error');
%%      pause();
%%  end

%%  %==============================
%%  % Angle Correction
%%  %==============================
%%  figure;
%%  surf(X, Y, grid.source(1).pixelAngleCorrection', 'EdgeColor', 'none');
%%  view(2);
%%  colorbar();
%%  figure;
%%  surf(X, Y, grid.source(10).pixelAngleCorrection', 'EdgeColor', 'none');
%%  view(2);
%%  colorbar();


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

