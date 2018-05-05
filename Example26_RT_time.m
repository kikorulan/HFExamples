% Compute the time propagation for the given grid
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex26_timeStep_imp;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex26_timeStep_imp;

close all;
%clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
load gridRT_Ray.mat;

%========================================
% Time propagation
%========================================
% Set time and forward signal
dt = grid.dt;
tMax = 2*dt;
grid.setTime(dt, tMax);
grid.setForwardSignal(1, [1 0 0]);

% Obtain 
grid.inverse_filter(100);
grid.timePropagationSignal(1);

%save gridRT_time.mat grid -v7.3;
%==================================================================================
% Plot results
%==================================================================================
x = 0:grid.dx:(grid.Nx-1)*grid.dx;
y = 0:grid.dy:(grid.Ny-1)*grid.dy;
  
%=========================
% Sound speed
%=========================
grid.plot_soundSpeed();
saveas(gcf, 'Example26_C.fig');
saveas(gcf, 'Example26_C', 'epsc');

%%  %=========================
%%  % Ray trajectories
%%  %=========================
%%  h = grid.plotRays(1, 10);

%%  %=========================
%%  % Distance surface
%%  %=========================
%%  figure;
%%  flipGrid = grid.source(1).pixelDistance';
%%  surf(x, y, flipGrid, 'EdgeColor', 'none');
%%  view(2);
%%  saveas(gcf, 'Example26_pixelDistance.fig');
%%  saveas(gcf, 'Example26_pixelDistance', 'epsc');
%%  
%%  %=========================
%%  % Time surface
%%  %=========================
%%  figure;
%%  flipGrid = grid.source(1).pixelTime';
%%  surf(x, y, flipGrid, 'EdgeColor', 'none');
%%  view(2);
%%  saveas(gcf, 'Example26_pixelTime.fig');

%=========================
% Attenuation surface
%=========================
figure;
flipGrid = grid.source(1).pixelAttenuation';
surf(x, y, log(abs(flipGrid)), 'EdgeColor', 'none');
view(2);
saveas(gcf, 'Example26_pixelAttenuation.fig');

%%  %=========================
%%  % Number of rays
%%  %=========================
%%  figure;
%%  flipGrid = pixelNRay';
%%  surf(x, y, flipGrid, 'EdgeColor', 'none');
%%  view(2);

%=========================
% Pressure surface
%=========================
meanRT = mean(abs(grid.source(1).pixelAPropagation(~isnan(grid.source(1).pixelAPropagation(:)))));
scrollView(permute(fliplr(real(grid.source(1).pixelAPropagation)), [2 1 3]), 3, [-6*meanRT 6*meanRT]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex26_timeStep_imp;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
