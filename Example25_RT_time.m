% Compute the time propagation for the given grid
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;

% Load filter
load Filter;
load gridRT_Ray;
%========================================
% Time propagation
%========================================

% Set time signal
sineSignal = sin(2*pi./tArray(end).*tArray);
grid.setTimeSignal(1, sineSignal, tArray);

%%% Time Varying Source Filter
grid.setCFilter(c0);
grid.setTFilter(tArray_TS);
grid.setFilter(Filter_TS);
grid_TS = grid.computeTimePropagation(1);
save grid_TS.mat grid_TS;

%%% Initial Value Source Filter
grid.setCFilter(c0);
grid.setTFilter(tArray_IV);
grid.setFilter(Filter_IV);
grid_IV = grid.computeTimePropagation(1);
save grid_IV.mat grid_IV nRays nSources;

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================

close all;
exampleFolder = '/cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/';
x = 0:grid.dx:(grid.Nx-1)*grid.dx;
y = 0:grid.dy:(grid.Ny-1)*grid.dy;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.c';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, strcat(exampleFolder, 'Ex25_timeStep/Example25_C.fig'));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%
h = grid.plotRays(1, 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance surface
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure;
flipGrid = grid.pixelDistance(:, :, 1)';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, strcat(exampleFolder, 'Ex25_timeStep/Example25_pixelDistance.fig'));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time surface
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.pixelTime(:, :, 1)';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, strcat(exampleFolder, 'Ex25_timeStep/Example25_pixelTime.fig'));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attenuation surface
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = grid.pixelAttenuation(:, :, 1)';
surf(x, y, log(abs(flipGrid)), 'EdgeColor', 'none');
view(2);
saveas(gcf, strcat(exampleFolder, 'Ex25_timeStep/Example25_pixelAttenuation.fig'));

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % Number of rays
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%  figure;
%%  flipGrid = pixelNRay';
%%  surf(x, y, flipGrid, 'EdgeColor', 'none');
%%  view(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure surface
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  meanRT = mean(abs(grid.pixelAmplitudeSignal(:)));
%%  scrollView(permute(fliplr(real(grid.pixelAmplitudeSignal)), [2 1 3]), 3, [-6*meanRT 6*meanRT])
