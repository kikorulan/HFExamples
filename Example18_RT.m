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
Nx = 120;
dx = 1e-3;
Ny = 180;
dy = 1e-3;
kernelSize = 33; % Multiple of 3

Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.2;
v2 = 0.2;
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

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];
Rgrid.setUMatrix(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays & sources
nRays = 2;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
x0 = [Nx*dx/2; Ny*dy/100]; % Source 1
% Create the new sources

% Sources

%%%% Proximal Method %%%
Rgrid.newSource(x0, 3*pi/7, 5*pi/7, nRays, step, tauMax);
Rgrid.computeSource(1);
Aprox = Rgrid.amplitudeProximal(1);
phiProx = Rgrid.source(1).phi;
qProx = Rgrid.source(1).q;

Rgrid.deleteSource(1);

%%%% ODE Method %%%
Rgrid.newSource(x0, 3*pi/7, 5*pi/7, nRays, step, tauMax);
Rgrid.computeSource(1);
Aode = Rgrid.revAmplitudeODE(1);
phiODE = Rgrid.source(1).phi;
qODE = Rgrid.source(1).q;

%%% New forward ray
[xEnd, angleEnd] = Rgrid.findTrajectoryEnd(1);
Rgrid.newSource(xEnd(:, :, 2), angleEnd(:, :, 2), angleEnd(:, :, 2), 1, step, tauMax);
Rgrid.computeSource(2);
Aforw = Rgrid.amplitudeProximal(2);
phiForw = Rgrid.source(2).phi;

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Speed of sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = Rgrid.c';
surf(x, y, flipGrid, 'EdgeColor', 'none');
colorbar;
view(2);
title('Speed of Sound');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex18_RT_interf_Q/Example18_C.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources and their respective Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
plot(Rgrid.source(1).x(1, :, 1), Rgrid.source(1).x(2, :, 1), 'Color', colourMapV(1), 'LineWidth', 2); % Main Ray 1
plot(Rgrid.source(1).x(1, :, 2), Rgrid.source(1).x(2, :, 2), 'Color', colourMapV(2), 'LineWidth', 2); % Main Ray 2
%plot(Rgrid.source(1).xPhi(1, :, 1), Rgrid.source(1).xPhi(2, :, 1), 'Color', colourMapV(1)); % Main Ray 1
%plot(Rgrid.source(1).xPhi(1, :, 2), Rgrid.source(1).xPhi(2, :, 2), 'Color', colourMapV(2)); % Main Ray 2
% Forward ray
plot(Rgrid.source(2).x(1, :, 1), Rgrid.source(2).x(2, :, 1), 'Color', colourMapV(3)); % Forward Ray
legend('Ray 1', 'Ray 2', 'Ray 2 - Forward');
title('Trajectories');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex18_RT_interf_Q/Example18_Ray_rev.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time vs Tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
plot(Rgrid.source(1).tau, Rgrid.source(1).phi(:, :, 1), 'Color', colourMapV(1)); % Main Ray 1
plot(Rgrid.source(1).tau, Rgrid.source(1).phi(:, :, 2), 'Color', colourMapV(2)); % Main Ray 2
legend('Ray 1', 'Ray 2');
title('Time');
xlabel('tau (m*m/s)');
ylabel('time (s)');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex18_RT_interf_Q/Example18_TimeVsTau_rev.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Proximal method
semilogy(phiProx(:, :, 1), Aprox(:, :, 1), 'Color', colourMapV(1), 'LineWidth', 2);  
hold on;
semilogy(phiProx(:, :, 2), Aprox(:, :, 2), 'Color', colourMapV(2), 'LineWidth', 2);  
% ODE method
semilogy(phiODE(:, :, 1), Aode(:, :, 1), 'Color', colourMapV(4));
semilogy(phiODE(:, :, 2), Aode(:, :, 2), 'Color', colourMapV(5));
% Forward ray
%semilogy(phiForw(:, :, 1), Aforw(:, :, 1), 'Color', colourMapV(3), 'LineWidth', 2);
legend('Ray 1 - Prox', 'Ray 2 - Prox', 'Ray 1 - ODE', 'Ray 2 - ODE', 'Ray 2 - Forward');
grid on;
xlabel('time (s)');
ylabel('Amplitude');
title('Amplitude (Reversal)');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex18_RT_interf_Q/Example18_A_rev.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Q - Determinant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proximal Method
figure;
hold on;
grid on;
plot(phiProx(:, :, 1), qProx(:, :, 1), 'Color', colourMapV(1));  
plot(phiProx(:, :, 2), qProx(:, :, 2), 'Color', colourMapV(2));  
legend('Ray 1', 'Ray 2');
xlabel('time (s)');
ylabel('q');
title('q (Prox reversal)');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex18_RT_interf_Q/Example18_Q_prox_rev.fig');

% ODE Method
figure;
hold on;
grid on;
plot(phiODE(:, :, 1), qODE(:, :, 1), 'Color', colourMapV(1));  
plot(phiODE(:, :, 2), qODE(:, :, 2), 'Color', colourMapV(2));  
legend('Ray 1', 'Ray 2');
xlabel('time (s)');
ylabel('q');
title('q (ODE reversal)');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex18_RT_interf_Q/Example18_Q_ODE_rev.fig');


