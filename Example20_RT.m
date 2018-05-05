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
Nx = 720;           % number of grid points in the x (row) direction
Ny = 1080;           % number of grid points in the y (column) direction
dx = 0.12/Nx;        % grid point spacing in the x direction [m]
dy = 0.18/Ny;        % grid point spacing in the y direction [m]
kernelSize = 33; % Multiple of 3

grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.1;
v2 = 0.1;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, floor(dimY/3));
M2 = v1*c0*ones(floor(dimX/2), floor(dimY/3));
M3 = v2*c0*ones(floor(dimX/2), floor(dimY/3));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [M1...
     [M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
grid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];
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
grid.newSource(x3, 0, 3*pi/8, nRays, step, tauMax);

% Sources
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    grid.computeSource(n);
    grid.computePressure(n);
end

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_C.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.u';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_U.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Sources
for n = 1:nSources
    % Rays
    for j = 1:nRays
        ray = grid.source(n).amplitude(:, :, j);
        tau = grid.source(n).tau;
        semilogy(tau(1:length(ray)), ray, 'Color', colourMapV(n));  
        hold on;
    end
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_Amplitude.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reverse Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Sources
for n = 1:nSources
    % Rays 
    for j = 1:nRays
        ray = grid.source(n).revAmplitude(:, :, j);
        tau = grid.source(n).tau;
        semilogy(tau(1:length(ray)), ray, 'Color', colourMapV(n));  
        hold on;
    end
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_RevAmplitude.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources and their respective Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
% Sources
for n = 1:nSources
    % Rays
    for j = 1:nRays
        plot(grid.source(n).x(1, :, j), grid.source(n).x(2, :, j), 'Color', colourMapV(n));
        %plot(grid.source(n).xPhi(1, :, j), grid.source(n).xPhi(2, :, j), 'Color', colourMapV(1));
    end
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_Ray.fig');

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
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_aSignal.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Filtered Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Spline
inc = 5e-8;
tSpline = 0:inc:1.3e-4;

maxFilter = find(Filter == max(Filter));
% Figure
figure;
hold on;
%axis([0 1.2e-4 -0.5 0.5]);
for n = 1:nSources
    lSpline = length(tSpline);
   
    tSignal = grid.source(n).tSignal;
    lSignal = length(tSignal(tSignal < inf));
    sSpline{n} = spline(grid.source(n).tSignal(1:lSignal), grid.source(n).aSignal(1:lSignal), tSpline);
    sConv = conv(sSpline{n}, Filter);
    plot(tSpline, sConv(1+maxFilter:lSpline+maxFilter)/max(Filter)*1e-3, 'Color', colourMapV(n));
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n+nSources));
end
legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex20_RT_kWave_interf/Example20_aSignalConv.fig');


end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
