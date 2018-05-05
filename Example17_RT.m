%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 480;
dx = 1e-3;
Ny = 720;
dy = 1e-3;
kernelSize = 33; % Multiple of 3

grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.25;
v2 = 0.25;
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
nRays = 2;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0*5/7;
step = min(dx, dy)*c0/12;

% Sources locations
x0 = [Nx*dx/2; Ny*dy/100]; % Source 1
% Create the new sources
grid.newSource(x0, 3*pi/7, 5*pi/7, nRays, step, tauMax);

% Sources
grid.computeSource(1);
grid.computeAmplitudeProximal(1);

%%%%
% Trace back ray
%%%%
% Sources
x1 = grid.source(1).x(:, end, 1);
x1C = grid.findCoordinates(x1);
c1 = grid.getC(x1C);
x2 = grid.source(1).x(:, end, 2);
% Angle
p1 = grid.source(1).p(:, end, 1);
p2 = grid.source(1).p(:, end, 2);
theta1 = atan(p1(2)/p1(1)) + pi*heaviside(-p1(1)) + pi;
theta2 = atan(p2(2)/p2(1)) + pi*heaviside(-p2(1)) + pi;
% Create the new sources
grid.newSource(x1, theta1, theta1, 1, step, tauMax);
grid.newSource(x2, theta2, theta2, 1, step, tauMax);

grid.computeSource(2);
grid.computeAmplitudeProximal(2);
grid.computeSource(3);
grid.computeAmplitudeProximal(3);


%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources and their respective Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
plot(grid.source(1).x(1, :, 1), grid.source(1).x(2, :, 1), 'Color', colourMapV(1)); % Main Ray 1
plot(grid.source(2).x(1, :, 1), grid.source(2).x(2, :, 1), 'Color', colourMapV(4)); % Back Ray 1
plot(grid.source(1).x(1, :, 2), grid.source(1).x(2, :, 2), 'Color', colourMapV(2)); % Main Ray 2
plot(grid.source(3).x(1, :, 1), grid.source(3).x(2, :, 1), 'Color', colourMapV(5)); % Back Ray 2
legend('Ray 1 - Forward', 'Ray 1 - Backward', ...
       'Ray 1 - Back Phi', 'Ray 2 - Forward', 'Ray 2 - Backward', ... 
       'MR 1', 'MR Phi 1', 'MR 2', 'MR Phi 2');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex17_RT_interf_Q/Example17_Ray.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Q & Amplitude - Source and End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q
figure;
hold on;
plot(grid.source(1).tau, grid.source(1).q(:, :, 1), 'Color', colourMapV(1));  
plot(grid.source(2).tau, grid.source(2).q(:, :, 1), 'Color', colourMapV(4));  
plot(grid.source(1).tau, grid.source(1).q(:, :, 2), 'Color', colourMapV(2));  
plot(grid.source(3).tau, grid.source(3).q(:, :, 1), 'Color', colourMapV(5));  
legend('Ray 1 - Forward', 'Ray 1 - Backward', 'Ray 2 - Forward', 'Ray 2 - Backward');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex17_RT_interf_Q/Example17_Q.fig');

% Amplitude
figure;
semilogy(grid.source(1).phi(:, :, 1), grid.source(1).amplitude(:, :, 1), 'Color', colourMapV(1));  
hold on;
semilogy(grid.source(2).phi(:, :, 1), grid.source(2).amplitude(:, :, 1), 'Color', colourMapV(4));  
semilogy(grid.source(1).phi(:, :, 2), grid.source(1).amplitude(:, :, 2), 'Color', colourMapV(2));  
semilogy(grid.source(3).phi(:, :, 1), grid.source(3).amplitude(:, :, 1), 'Color', colourMapV(5));  
legend('Ray 1 - Forward', 'Ray 1 - Backward', 'Ray 2 - Forward', 'Ray 2 - Backward');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex17_RT_interf_Q/Example17_A.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Interface loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I1 = grid.interfaceLoss(1);
disp('Hola');
I2 = grid.interfaceLoss(2);
I3 = grid.interfaceLoss(3);

figure;
hold on;
plot(grid.source(1).tau, I1(:, :, 1), 'Color', colourMapV(1));
plot(grid.source(2).tau, fliplr(I2(:, :, 1)), 'Color', colourMapV(4));
plot(grid.source(1).tau, I1(:, :, 2), 'Color', colourMapV(2));
plot(grid.source(3).tau, fliplr(I3(:, :, 1)), 'Color', colourMapV(5));
legend('Ray 1 - Forward', 'Ray 1 - Backward', 'Ray 2 - Forward', 'Ray 2 - Backward');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex17_RT_interf_Q/Example17_Interface.fig');
