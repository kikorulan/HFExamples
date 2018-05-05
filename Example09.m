%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Measure computational time
tic;
start_time = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 120;
dx = 1e-3;
Ny = 120;
dy = 1e-3;
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
% Auxiliary matrices
c = c0*ones(Nx, Ny);
grid.setCMatrix(c);

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

% Number of rays
nRays = 10;

tauMax = 2*(Nx*dx+Ny*dy)*c0;
step = min(dx, dy)*c0/3;

%%%%%%
amplitude = 1;
x0 = [[Nx*dx/2]; [Ny*dy/10]]; % Ray 1
% Run the phase
%grid.computePhase(x0);

% Create a new source
grid.newSource(x0, step, tauMax);

thetaMin = pi/3;
thetaMax = 2*pi/3;
% Rays
for j = 1:nRays
    theta = thetaMin + (j-1)*(thetaMax-thetaMin)/nRays;
    p0 = [cos(theta); sin(theta)]/c0;
    grid.computeRay(1, p0);
    % Reverse Amplitude
    grid.computeReverseAmplitude(1, j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 1:1:Nx;
y = 1:1:Ny;

% Rays
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
for i = 1:1:nRays
    plot(grid.source(1).ray(i).x(1, :), grid.source(1).ray(i).x(2, :));
    plot(grid.source(1).ray(i).xPhi(1, :), grid.source(1).ray(i).xPhi(2, :));
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex09_RT_ampl/Example9_Ray.fig');

% N
figure;
flipGrid = grid.u';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex09_RT_ampl/Example9_U.fig');

% Amplitude
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy 0 3]);
for j = 1:nRays
    lengthRay = length(grid.source(1).ray(j).revAmplitude);
    ray = {grid.source(1).ray(j).x(1, 1:lengthRay), grid.source(1).ray(j).x(2, 1:lengthRay), grid.source(1).ray(j).revAmplitude};
    plot3(ray{1}, ray{2}, ray{3});  
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex09_RT_ampl/Example9_RevAmplitude.fig');

% Amplitude
figure;
lengthA = length(grid.source(1).ray(1).amplitude);
plot(grid.source(1).ray(1).tau(1:lengthA), grid.source(1).ray(1).amplitude);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex09_RT_ampl/Example9_Amplitude.fig');


end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

