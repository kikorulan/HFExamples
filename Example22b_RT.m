%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;

close all;
%clear all;

load kgrid_data.mat;

% Measure computational time
tic;
start_time = clock;
 
run colourMap;
%========================================
% Grid definition
%========================================
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]

grid = gridRT(Nx, dx, Ny, dy);
% Sound speed
grid.setCMatrix(cMatrix);
% Initial pressure
grid.setUMatrix(source.p0);
%========================================
% Impulse Response
%========================================
% Set time
dt = min(grid.dx, grid.dy)/c0/3;
tMax = 4e-5;
grid.setTime(dt, tMax);

grid.impulse_additive('IV');

save gridRT_impulse.mat grid;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse.mat;
% Number of rays & sources
nRays = 1000;
nSources = 3;

tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*c0;
tauStep = min(grid.dx, grid.dy)*c0/3;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/c0;
tStep = min(grid.dx, grid.dy)/c0/3;

% Sources locations
x1(1, 1, 1) = grid.Nx*grid.dx/2; 
x1(1, 1, 2) = grid.Ny*grid.dy/100; 
x2(1, 1, 1) = 5*grid.Nx*grid.dx/6; 
x2(1, 1, 2) = grid.Ny*grid.dy/100; 
x3(1, 1, 1) = grid.Nx*grid.dx/100; 
x3(1, 1, 2) = grid.Ny*grid.dy/2; 


% Create the new sources
grid.setDeltaX(grid.dx/100);
clear source;
source(1) = grid.newSource(x1, 0, pi, nRays, tStep, tMax);
source(2) = grid.newSource(x2, 0, pi, nRays, tStep, tMax);
source(3) = grid.newSource(x3, -pi/2, pi/2, nRays, tStep, tMax);
  
% Sources
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    grid.computeHamil(source(n));
end

% Save results
%save gridRT_ray.mat grid nRays nSources x1 x2 x3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;

%load gridRT.mat;
%load kgrid_data.mat;

x = 0:grid.dx:(grid.Nx-1)*grid.dx;
y = 0:grid.dy:(grid.Ny-1)*grid.dy;

%%  %==============================
%%  % Sound speed
%%  %==============================
%%  grid.plotSoundSpeed();
%%  saveas(gcf, 'Example22b_C', 'epsc'); 
%%  saveas(gcf, 'Example22b_C.fig'); 

%%  %==============================
%%  % Initial pressure
%%  %==============================
%%  figure;
%%  flipGrid = grid.u';
%%  surf(x, y, flipGrid, 'EdgeColor', 'none');
%%  view(2);
%%  saveas(gcf, 'Example22b_U', 'epsc');
%%  saveas(gcf, 'Example22b_U.fig'); 

%==============================
% Sources and their respective Rays
%==============================
figure;
hold on;
axis([0 grid.Nx*grid.dx 0 grid.Ny*grid.dy]);
% Sources
for n = 1:nSources
    % Rays
    for j = 1:10:nRays
        plot(source(n).x(j, :, 1), source(n).x(j, :, 2), 'Color', colourMapV(n));
    end
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, 'Example22b_Ray', 'epsc');
saveas(gcf, 'Example22b_Ray.fig');

%==============================
% Non-Filtered Signal
%==============================
% Time Signal Amplitude
figure;
hold on;
for n = 1:nSources
    plot(source(n).tBeam, source(n).aBeam, 'Color', colourMapV(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3')
saveas(gcf, 'Example22b_aSignal', 'epsc');
saveas(gcf, 'Example22b_aSignal.fig');

%==============================
% Filtered Signal
%==============================
% Normalisation factor
normRT = max(source(1).aForward);
normKWave = max(sensor_data.p(1, :));

% Figure
figure;
hold on;
%axis([0 1.2e-4 -0.5 0.5]);
for n = 1:nSources
    % Plot and compare
    plot(grid.tForward, source(n).aForward/normRT, 'Color', colourMapV(n));
    plot(kgrid.t_array, sensor_data.p(n, :)/normKWave, 'Color', colourMapV(n+nSources), 'LineWidth', 2);
end

legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
saveas(gcf, 'Example22b_aSignalConv', 'epsc');
saveas(gcf, 'Example22b_aSignalConv.fig');

%==============================
% Isotime curve
%==============================
nCurves = 20;
nS = 3;
xIso = grid.findIsoTime(nS, 20, 0.6e-5, 1.4e-5);
% Figure
figure;
flipGrid = grid.u';
surf(x, y, flipGrid, 'EdgeColor', 'none');
hold on;
for n = 1:nCurves
    plot3(xIso(n, :, 1), xIso(n, :, 2), repmat(2, [1 size(xIso, 2)]), '-m');
end
for j = 1:20:nRays
    plot3(source(nS).x(j, :, 1), source(nS).x(j, :, 2), repmat(2, [1 size(grid.source(nS).x, 2)]), '-r');
end
view(2);
legend('Isocurves for times tMin = 0.9e-5s, tMax = 1.1e-5s');
saveas(gcf, 'Example22b_isotime.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/;


