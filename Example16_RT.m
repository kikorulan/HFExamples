%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

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

Rgrid = gridRT(Nx, dx, Ny, dy);

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

Rgrid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of sources
nRays = 5;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/12;

% Source
x0 = [[Nx*dx/2]; 20*dy];

% Create the new source
Rgrid.newSource(x0, pi/2+0.4, pi/2, nRays, step, tauMax);
Rgrid.computeSource(1);
A = Rgrid.computeAmplitudeSI(1);
B = Rgrid.computeAmplitudeProximal(1);


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
    xCoord = Rgrid.source(1).x(1, :, i);
    yCoord = Rgrid.source(1).x(2, :, i);
    plot(xCoord, yCoord, 'Color', colourMapBlue(i));
end
title('Rays');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex16_RT_interf/Example16_Rays.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = Rgrid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex16_RT_interf/Example16_C.fig'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
axis([0 7e-4 1e-4 1]);
% Rays
for j = 1:nRays
    phi = Rgrid.source(1).phi(:, :, j);
    semilogy(phi(1:length(A(:, :, j))), A(:, :, j), 'Color', colourMapRed(j), 'LineWidth' , 2);
    hold on;
end
for j = 1:nRays
    phi = Rgrid.source(1).phi(:, :, j);
    semilogy(phi(1:length(B(:, :, j))), B(:, :, j), 'Color', colourMapBlue(j), 'LineWidth' , 1);
end
title('Amplitude');
xlabel('t');
ylabel('Amplitude');
legend('Ray 1 - SI', 'Ray 2 - SI', 'Ray 3 - SI', 'Ray 4 - SI', 'Ray 5 - SI', ...
       'Ray 1 - PR', 'Ray 2 - PR', 'Ray 3 - PR', 'Ray 4 - PR', 'Ray 5 - PR');
grid on;
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex16_RT_interf/Example16_Amplitude.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Traveled distance
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Rays
for j = 1:nRays
    phi = Rgrid.source(1).phi(:, :, j);
    plot(phi, Rgrid.source(1).xTD(:, :, j), 'Color', colourMapRed(j), 'LineWidth' , 2);
    hold on;
end
title('Traveled distance');
xlabel('t');
ylabel('Traveled Distance');
grid on;
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex16_RT_interf/Example16_TD.fig'); 

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
