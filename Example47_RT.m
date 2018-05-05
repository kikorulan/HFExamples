%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex47_GB;

close all;
%clear all;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

% Make matrix symmetric
c = 1500*ones(Nx, Ny);
gridR.setCMatrix(c);
% Set initial pressure
u = zeros(Nx, Ny);
gridR.setUMatrix(u);

%%  %========================================
%%  % Impulse Response
%%  %========================================
%%  % Set time
%%  dt = 2e-7;
%%  %dt = min(gridR.dx, gridR.dy)/c0/2;
%%  tMax = 2e-4;
%%  gridR.setTime(dt, tMax);
%%  % Compute impulse response
%%  gridR.impulse_additive('IV');
%%  save gridRT_impulse.mat gridR dt tMax;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 1;% 800
nSources = 1;%256

% Parametrisation
tStep = dt;

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-1)/2*gridR.dx, 0);

% Sources
source(1) = gridR.newSource(x{1}, pi/4, 3*pi/4, nRays, tStep, tMax);
gridR.computeHamil(source(1));

%  clear source;
%%  source = gridR.computeForwardParallel(x, 0, pi, nRays, tStep, tMax, false);

% Save results
%save gridRT.mat gridR nRays nSources x -v7.3;
  
%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex47_GB;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');


%==============================
% Amplitude
%==============================
figure;
plot(gridR.tForward, real(source(1).amplitude), 'Color', 'r');
hold on;
plot(gridR.tForward, imag(source(1).amplitude), 'Color', 'b');
legend('amplitude - real', 'amplitude - imag');
title('Amplitudes');
xlabel('t (s)');

%saveas(gcf, 'Example46_amplitudeODE_log_caustic.fig');

%==============================
% Q & P
%==============================

% Q
figure;
plot(gridR.tForward, real(source(1).q), 'Color', 'r');
hold on;
plot(gridR.tForward, imag(source(1).q), 'Color', 'b');
legend('q - real', 'q - imag');
grid on;
box on;
title('Qs');
xlabel('t (s)');

% P
figure;
plot(gridR.tForward, real(source(1).kIndex), 'Color', 'r');
hold on;
plot(gridR.tForward, imag(source(1).kIndex), 'Color', 'b');
legend('p - real', 'p - imag');
grid on;
box on;
title('Ps');
xlabel('t (s)');

%==============================
% Curvature and Beamwidth
%==============================

% Curvature
curvature = 1./source(1).n.*real(source(1).kIndex./source(1).q);
figure;
plot(gridR.tForward, curvature(1, :), 'Color', 'r');
grid on;
box on;
title('Curvature');
xlabel('t (s)');

% Beamwidth
beamwidth = sqrt(-1./imag(source(1).kIndex./source(1).q));
figure;
plot(gridR.tForward, beamwidth, 'Color', 'r');
grid on;
box on;
title('beamwidth');
xlabel('t (s)');


%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

