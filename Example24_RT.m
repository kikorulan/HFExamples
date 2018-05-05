%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

load('workspace');
% Measure computational time
tic;
start_time = clock;

run colourMap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 200;           % number of grid points in the x (row) direction
Ny = 200;           % number of grid points in the y (column) direction
dx = 0.1/Nx;        % grid point spacing in the x direction [m]
dy = 0.1/Ny;        % grid point spacing in the y direction [m]

grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;

% Randomise C
c = c0*ones(Nx, Ny);
grid.setCMatrix(c);
spatialSDV = 0.09;
deltaC = 0.1;
grid.randomiseC(spatialSDV, deltaC);

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/4),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M2 ; M3] M4];
grid.setUMatrix(U);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of rays & sources
nRays = 300;
nSources = 3;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/6;

% Sources locations
x1 = [Nx*dx/2; Ny*dy/100]; % Source 1
x2 = [5*Nx*dx/6; Ny*dy/100]; % Source 2
x3 = [Nx*dx/100; Ny*dy/3]; % Source 3

% Create the new sources
grid.newSource(x1, pi/3, 2*pi/3, nRays, step, tauMax);
grid.newSource(x2, pi/2, 5*pi/6, nRays, step, tauMax);
grid.newSource(x3, 0, 15*pi/16, nRays, step, tauMax);

% Sources
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    grid.computeSourceQ(n);
    %A = grid.amplitudeODE(n);
    grid.computePressure(n);
end

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
close all;
x = 1:1:Nx;
y = 1:1:Ny;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.c';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_randomC/Example24_C.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.u';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_randomC/Example24_U.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 grid.Nx*grid.dx 0 grid.Ny*grid.dx]);
for n = 1:nSources
    for j = 1:nRays
        plot(grid.source(n).x(1, :, j), grid.source(n).x(2, :, j), 'Color', colourMapV(n));
    end
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_randomC/Example24_Rays.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%
h = grid.plotAmplitude(1, 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Non-Filtered Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Signal Amplitude
figure;
hold on;
for n = 1:nSources
    grid.timeSignal(n);
    plot(grid.source(n).tSignal, grid.source(n).aSignal, 'Color', colourMapV(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_randomC/Example24_aSignal.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Filtered Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Spline
tSpline = 0:dt:max(kgrid.t_array);
lSpline = length(tSpline);

% Figure
figure;
hold on;
for n = 1:nSources
    maxFilter{n} = find(Filter{n} == max(Filter{n}));
    tSignal = grid.source(n).tSignal;
    lSignal = length(tSignal(tSignal < inf));
    sSpline{n} = spline(grid.source(n).tSignal(1:lSignal), grid.source(n).aSignal(1:lSignal), tSpline);
    sConv = conv(sSpline{n}, Filter{n});
    plot(tSpline, sConv(1+maxFilter{n}:lSpline+maxFilter{n})/max(Filter{n})*0.2, 'Color', colourMapV(n));
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n+nSources));
end
legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_randomC/Example24_aSignalConv.fig');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

