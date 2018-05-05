%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex44_WENO;

close all;
clear all;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 200;         % number of grid points in the x (row) direction
Ny = 200;         % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
grid.setCMatrix(c0*ones(Nx, Ny));

% Set initial pressure
grid.setUMatrix(zeros(Nx, Ny));


%%  %========================================
%%  % Impulse Response
%%  %========================================
%%  % Set time
%%  dt = 3e-7;
%%  %dt = min(grid.dx, grid.dy)/c0/2;
%%  tMax = 2e-4;
%%  grid.setTime(dt, tMax);
%%  % Compute impulse response
%%  grid.impulse_additive('IV');
%%  
%%  save gridRT_impulse.mat grid c0 dt;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse.mat;
% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 1000;% 800
nSources = 1;

% Parametrisation
tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/c0;
tStep = dt/10;

% Sources locations
clear x;
x{1} = cat(3, (Nx-1)/2*dx, 0);

clear source;
source = grid.computeForwardParallel(x, 0, pi, nRays, tStep, tMax, false);

% Save aForward data
% Assign data
aForward_RT = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    aForward_RT(n, :) = source(n).aForward;
end

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex44_WENO;

%load gridRT.mat;
%load sensor_data.mat;


%==============================
% Sources and their respective Rays
%==============================
figure;
axis([0 grid.xAxis(end) 0 grid.yAxis(end)]);
hold on;
colours = winter(nRays);
for n = 1:nRays
    plot(source.x(n, :, 1), source.x(n, :, 2), 'Color', colours(n, :));
end

%==============================
% Propagation time
%==============================
pixelTime = source(1).pixelTime';
figure;
contour(pixelTime, 30);

%Derivative with respect to x
figure;
pixelTime_dx = pixelTime(2:end, :)-pixelTime(1:end-1, :);
surf(grid.xAxis, grid.yAxis(1:end-1), pixelTime_dx, 'EdgeColor', 'none');
view(2);
figure;
contour(grid.xAxis, grid.yAxis(1:end-1), pixelTime_dx, 30);
xlabel('x (m)');
ylabel('y (m)');
title('Derivative of time with respect to x (coarse mesh)');
saveas(gcf, 'Example44_derivative_coarse.fig');

%==============================
% Amplitude attenuation
%==============================
pixelAttenuation = source(1).pixelAttenuation';
figure;
contour(grid.xAxis, grid.yAxis, log(pixelAttenuation), 30);
xlabel('x (m)');
ylabel('y (m)');
title('Amplitude attenuation (coarse mesh) - log scale');
saveas(gcf, 'Example44_amplitude_coarse.fig');

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

