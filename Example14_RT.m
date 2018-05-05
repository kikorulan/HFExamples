%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% Measure computational time
tic;
start_time = clock;

run colorMap;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of sources
nRays = 40;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Source
x0 = [[Nx*dx/3]; 20*dy];

% Create the new source
grid.newSource(x0, pi/2, pi/2+0.4, nRays, step, tauMax);
grid.computeSource(1);
grid.computeAmplitude_ODE(1);

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
for i = 1:1:nRays
    xCoord = grid.source(1).x(1, :, i);
    yCoord = grid.source(1).x(2, :, i);
    plot(xCoord, yCoord, 'Color', colourMapBlue(i));
end
title('Rays');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex14_RT_interf/Example14_Rays.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);

saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex14_RT_interf/Example14_C.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation times
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
% Rays
tau = grid.source(1).tau;
for j = 1:nRays
    phi = grid.source(1).phi(:, :, j);
    plot(tau, phi, 'Color', colourMapBlue(j));
end
title('Tau');
xlabel('tau');
ylabel('time');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex14_RT_interf/Example14_Time.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude - Delta X
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
semilogy(0, 1);
axis([0 7e-4 1e-4 1]);
hold on;
% Rays
for j = 1:nRays
    ray = grid.source(1).amplitude(:, :, j);
    phi = grid.source(1).phi(:, :, j);
    semilogy(phi(1:length(ray)), ray, 'Color', colourMapBlue(j));
end
title('Amplitude');
xlabel('t');
ylabel('Amplitude');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex14_RT_interf/Example14_Amplitude.fig'); 

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
