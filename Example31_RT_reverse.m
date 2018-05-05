%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

close all;
%clear all;

%load gridRT.mat;
load recon_data_sensors.mat;
load recon_data_sensors_TVsource.mat;
% Measure computational time
tic;
start_time = clock;
 
%run colourMap;
 
%========================================
% Compute filters
%========================================
grid.inverse_filter(10);
%grid.impulse_dirichlet();
%========================================
% Compute reverse signal
%========================================
%grid.inverse_beam_adjoint();

for n = 1:nSources
    grid.inverse_beam(source(n));
end

%grid.computeTimeReverseParallel(source);
pixelAReverseSensors = grid.inverse_signal(source);


%save gridRT_recon.mat grid -v7.3;
%================================================================================
% VISUALISATION
%================================================================================
X = 0:grid.dx:(grid.Nx-1)*grid.dx;
Y = 0:grid.dy:(grid.Ny-1)*grid.dy;
axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 400 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%%  %========================================
%%  % Forward Signal
%%  %========================================
%%  figure;
%%  plot(grid.tForward, grid.source(1).aForward);
%%  title('Forward Signal');

%%  %========================================
%%  % Reverse Filter
%%  %========================================
%%  figure;
%%  plot(grid.tFilter, grid.FilterReverse);
%%  title('Reverse Filter');

%%  %========================================
%%  % Reverse Signal
%%  %========================================
%%  figure;
%%  plot(grid.tReverse, grid.source(1).aReverse);
%%  title('Reverse Signal');

%==============================
% Initial Pressure
%==============================
h = figure;
surf(X, Y, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Initial Pressure');
saveas(gcf, 'Example31_U','png'); 
saveas(gcf, 'Example31_U.fig'); 

%========================================
% Reconstruction - RT
%========================================
figure;
surf(X, Y, grid.pixelAReverse', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction - Ray Tracing');
saveas(gcf, 'Example31_RT_recon_ADJ', 'png'); 
saveas(gcf, 'Example31_RT_recon_ADJ.fig');

%========================================
% Reconstruction - kWave (ADJ)
%========================================
%pixelKWave = p0_recon{257}/max(abs(p0_recon{257}(:)));
pixelKWave_ADJ = max(0, p0_recon_TVsource{257}.p_final);
pixelKWave_ADJ = pixelKWave_ADJ/max(pixelKWave_ADJ(:));
figure;
surf(X, Y, pixelKWave_ADJ', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction - kWave (ADJ)');
saveas(gcf, 'Example31_kWave_recon_ADJ', 'png'); 
saveas(gcf, 'Example31_kWave_recon_ADJ.fig');

%========================================
% Error - RT (ADJ)
%========================================
errorRT_ADJ = pixelKWave_ADJ - grid.pixelAReverse;
figure;
surf(X, Y, errorRT_ADJ', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Error - Ray Tracing (ADJ)');
saveas(gcf, 'Example31_RT_error_ADJ', 'png'); 
saveas(gcf, 'Example31_RT_error_ADJ.fig');

%========================================
% Reconstruction - kWave 
%========================================
%pixelKWave = p0_recon{257}/max(abs(p0_recon{257}(:)));
pixelKWave_TR = max(0, p0_recon{257});
pixelKWave_TR = pixelKWave_TR/max(pixelKWave_TR(:));
figure;
surf(X, Y, pixelKWave_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction - kWave (TR)');
saveas(gcf, 'Example31_kWave_recon_TR', 'png'); 
saveas(gcf, 'Example31_kWave_recon_TR.fig');

%========================================
% Error - RT
%========================================
errorRT_TR = pixelKWave_TR - grid.pixelAReverse;
figure;
surf(X, Y, errorRT_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Error - Ray Tracing (TR)');
saveas(gcf, 'Example31_RT_error_TR', 'png'); 
saveas(gcf, 'Example31_RT_error_TR.fig');

%%  %========================================
%%  % Error - kWave
%%  %========================================
%%  errorKWave = grid.u - pixelKWave;
%%  figure;
%%  surf(X, Y, errorKWave', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  caxis([-1, 1]);
%%  colorbar();
%%  title('Error - kWave');
%%  saveas(gcf, 'Example31_kWave_error.fig');
%%  
%%  %========================================
%%  % Reconstruction - Sensor 1
%%  %========================================
%%  figure;
%%  pixelAReverse = grid.source(128).pixelAReverse;
%%  pixelAReverse(isnan(pixelAReverse)) = 0;
%%  pixelAReverse = pixelAReverse/max(pixelAReverse(:));
%%  surf(X, Y, pixelAReverse', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  colorbar();
%%  box on;
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  set(gcf, 'pos', position);
%%  %title('Reconstruction - Sensor 128');
%%  saveas(gcf, 'Example31_RT_sensor128', 'png'); 
%%  saveas(gcf, 'Example31_RT_sensor128.fig');

%%  %========================================
%%  % Reconstruction - Sensors
%%  %========================================
%%  meanSensors = mean(abs(pixelAReverseSensors(:)));
%%  scrollView(permute(fliplr(real(pixelAReverseSensors)), [2 1 3]), 3, [-6*meanSensors 6*meanSensors])

%%  %========================================
%%  % Pixel Time propagation
%%  %========================================
%%  figure;
%%  surf(X, Y, grid.pixelTime', 'EdgeColor', 'none');
%%  view(2);
%%  colorbar();
%%  title('Time Propagation');

%%  %========================================
%%  % Pixel Reverse Time
%%  %========================================
%%  figure;
%%  surf(X, Y, grid.pixelTReverse', 'EdgeColor', 'none');
%%  view(2);
%%  colorbar();
%%  title('Reverse Time');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
