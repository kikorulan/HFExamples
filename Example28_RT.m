% Example for gridRT class
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex28_Shadow;

close all;
clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;

%========================================
% Grid definition
%========================================
Nx = 720;           % number of grid points in the x (row) direction
Ny = 720;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
kernelSize = 33; % Multiple of 3

Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -.2;
v2 = .2;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, floor(dimY/2));
M2 = v1*c0*ones(floor(dimX/2), floor(dimY/2));
M3 = v2*c0*ones(floor(dimX/2), floor(dimY/2));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [[M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
Rgrid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

%========================================
% Impulse Response
%========================================
% Set time
dt = 3e-7;
tMax = 6e-4;
Rgrid.setTime(dt, tMax);
Rgrid.impulseResponse2D('IV');

save gridRT_impulse.mat Rgrid;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse.mat;

% Number of rays & sources
nRays1 = 50;
nRays2 = 200;

tMax = 8e-4;
tStep = min(Rgrid.dx, Rgrid.dy)/c0/3;

% Sources locations
x1(1, 1, 1) = Rgrid.Nx*Rgrid.dx*(0.3); 
x1(1, 1, 2) = (Rgrid.Ny-1)*Rgrid.dy; 
x2(1, 1, 1) = 0;
x2(1, 1, 2) = Rgrid.Ny*Rgrid.dy*0.1;
% Create the new sources
Rgrid.newSource(x1, -pi/2, -pi/4, nRays1, tStep, tMax);
Rgrid.newSource(x2, pi/32, pi/16, nRays2, tStep, tMax);

% Sources
Rgrid.computeHamil(1);
Rgrid.computeHamil(2);

save gridRT.mat Rgrid;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex28_Shadow;

load gridRT;
%========================================
% Sound Speed
%========================================
Rgrid.plotSoundSpeed();
title('');
box on;
set(gcf, 'pos', [100 100 600 500])
saveas(gcf, 'Example28_C', 'png');
saveas(gcf, 'Example28_C.fig');

%========================================
% Caustics
%========================================
% Ray Trajectories
Rgrid.plotRays(1, 100);
C1 = plot(0.36, 0.48, 'Color', 'r', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 3, 'LineStyle', 'none');
C2 = plot(0.35, 0.40, 'Color', 'g', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 3, 'LineStyle', 'none');
C3 = plot(0.28, 0.29, 'Color', 'r', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 3, 'LineStyle', 'none');
C4 = plot(0.20, 0.22, 'Color', 'g', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 3, 'LineStyle', 'none');
C5 = plot(0.14, 0.16, 'Color', 'g', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 3, 'LineStyle', 'none');
C6 = plot(0.08, 0.10, 'Color', 'g', 'Marker', 'o', 'MarkerSize', 20, 'LineWidth', 3, 'LineStyle', 'none');
leg = legend([C1 C2], 'Caustic type 1', 'Caustic type 2');
HeightScaleFactor = 1.7;
NewHeight = leg.Position(4) * HeightScaleFactor;
leg.Position(2) = leg.Position(2) - (NewHeight - leg.Position(4));
leg.Position(4) = NewHeight;
grid on;
box on;
title('');
set(gcf, 'pos', [700 700 500 500]);
saveas(gcf, 'Example28_Caustics', 'epsc');
saveas(gcf, 'Example28_Caustics.fig');

% Amplitude
Rgrid.plotAmplitude(1, 100);
title('');
saveas(gcf, 'Example28_amplitude', 'epsc');
saveas(gcf, 'Example28_amplitude.fig');
% Q
Rgrid.plotQ(1, 100);
axis([0 8e-4 -4e-6 5e-6]);
title('');
saveas(gcf, 'Example28_Q', 'epsc');
saveas(gcf, 'Example28_Q.fig');

%========================================
% Caustics
%========================================
% Ray Trajectories
Rgrid.plotRays(2, 100);
plot(0.48, 0.6, 'Color', 'r', 'Marker', 's', 'MarkerSize', 70, 'LineWidth', 3);
grid on;
box on;
title('');
set(gcf, 'pos', [700 700 500 500]);
saveas(gcf, 'Example28_ShadowRegions', 'epsc');
saveas(gcf, 'Example28_ShadowRegions.fig');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

