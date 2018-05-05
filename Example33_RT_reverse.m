%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex33_reconstruction;

close all;
%clear all;

%load gridRT.mat;
%load recon_data.mat;
load recon_data_adjoint.mat;
% Measure computational time
tic;
start_time = clock;
 
%run colourMap;
 
%========================================
% Compute filters
%========================================
grid.inverse_filter(100);

%========================================
% Compute reverse signal
%========================================
%grid.inverse_beam(1);

%aReverse = grid.inverse_beam_adjoint();
for n = 1:nSources
    grid.inverse_beam(source(n));
end
%grid.computeTimeReverseParallel();

pixelAReverseSensors = grid.inverse_signal(source);

%========================================
% Subsample versions of the reconstruction
%========================================
pixelAReverse = grid.pixelAReverse;
% 128 sensors
source_sub_128 = source(1:2:end);
grid.inverse_signal(source_sub_128);
pixelAReverse_128 = grid.pixelAReverse;

% 64 sensors
source_sub_64 = source(2:4:end);
grid.inverse_signal(source_sub_64);
pixelAReverse_64 = grid.pixelAReverse;

% 32 sensors
source_sub_32 = source(4:8:end);
grid.inverse_signal(source_sub_32);
pixelAReverse_32 = grid.pixelAReverse;

% 16 sensors
source_sub_16 = source(8:16:end);
grid.inverse_signal(source_sub_16);
pixelAReverse_16 = grid.pixelAReverse;

save recon_data_RT.mat pixelAReverse pixelAReverse_128 pixelAReverse_64 pixelAReverse_32 pixelAReverse_16;
%================================================================================
% VISUALISATION
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex33_reconstruction;

axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 320 600];
positionNoY = [700 700 300 600];
positionNoYBar = [700 700 363 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% Sound Speed
%========================================
grid.plot_soundSpeed();
title('');
box on;
set(gcf, 'pos', positionYBar);
%title('Sound Speed');
saveas(gcf, 'Example33_C', 'png'); 
saveas(gcf, 'Example33_C.fig'); 

%==============================
% Initial Pressure
%==============================
figure;
surf(grid.xAxis, grid.yAxis, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionYBar);
%title('Initial Pressure');
saveas(gcf, 'Example33_U', 'png'); 
saveas(gcf, 'Example33_U.fig'); 

%========================================
% Reconstruction - RT
%========================================
maxPixelRT = max(real(pixelAReverse(:)));
pixelRT = max(0, real(pixelAReverse)/maxPixelRT);
figure;
surf(grid.xAxis, grid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction RT');
saveas(gcf, 'Example33_RT_recon', 'png');
saveas(gcf, 'Example33_RT_recon.fig');

%========================================
% Reconstruction - kWave Adjoint
%========================================
pixelKWave = max(0, p0_recon_adjoint.p_final/max(p0_recon_adjoint.p_final(:)));
figure;
surf(grid.xAxis, grid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
%title('Reconstruction k-Wave');
saveas(gcf, 'Example33_kWave_recon', 'png');
saveas(gcf, 'Example33_kWave_recon.fig');

%========================================
% Error - kWave-RT
%========================================
errorRT = pixelKWave - pixelRT;
figure;
surf(grid.xAxis, grid.yAxis, errorRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
%title('Error RT - kWave (ADJ)');
saveas(gcf, 'Example33_RT_error', 'png');
saveas(gcf, 'Example33_RT_error.fig');

%%  %========================================
%%  % Error - kWave
%%  %========================================
%%  errorKWave = grid.u - pixelKWave;
%%  figure;
%%  surf(grid.xAxis, grid.yAxis, errorKWave', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  caxis([-1, 1]);
%%  colorbar();
%%  title('Error: kWave');
%%  saveas(gcf, 'Example33_kWave_error.fig');

%%  %========================================
%%  % Reconstruction - Sensor 1
%%  %========================================
%%  figure;
%%  pixelAReverse = real(source(128).pixelAReverse);
%%  pixelAReverse(isnan(pixelAReverse)) = 0;
%%  pixelAReverse = pixelAReverse/max(pixelAReverse(:));
%%  surf(grid.xAxis, grid.yAxis, pixelAReverse', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  colorbar();
%%  box on;
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  set(gcf, 'pos', position);
%%  %title('Reconstruction - Sensor 128');
%%  saveas(gcf, 'Example33_RT_sensor128', 'png'); 
%%  saveas(gcf, 'Example33_RT_sensor128.fig');

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

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
