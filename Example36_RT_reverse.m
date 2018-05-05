%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;

close all;
%clear all;

%load gridRT.mat;
load recon_data.mat;
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



%save gridRT_recon.mat grid -v7.3;
%================================================================================
% VISUALISATION
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;

axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 400 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% Sound Speed
%========================================
grid.plot_soundSpeed();
title('');
box on;
set(gcf, 'pos', position);
title('Sound Speed');
saveas(gcf, 'Example36_C', 'png'); 
saveas(gcf, 'Example36_C.fig'); 

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
figure;
surf(grid.xAxis, grid.yAxis, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Initial Pressure');
saveas(gcf, 'Example36_U', 'png'); 
saveas(gcf, 'Example36_U.fig'); 

%========================================
% Reconstruction - RT
%========================================
maxPixelRT = max(real(grid.pixelAReverse(:)));
pixelRT = max(0, real(grid.pixelAReverse)/maxPixelRT);
figure;
surf(grid.xAxis, grid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction RT');
saveas(gcf, 'Example36_RT_recon', 'png');
saveas(gcf, 'Example36_RT_recon.fig');

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
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction k-Wave (ADJ)');
saveas(gcf, 'Example36_kWave_recon', 'png');
saveas(gcf, 'Example36_kWave_recon.fig');

%========================================
% Error - RT-kWave
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
ylabel('y (m)');
set(gcf, 'pos', position);
title('Error RT - kWave (ADJ)');
saveas(gcf, 'Example36_RT_error', 'png');
saveas(gcf, 'Example36_RT_error.fig');

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
%%  saveas(gcf, 'Example36_kWave_error.fig');

%========================================
% Reconstruction - Sensor 1
%========================================
figure;
pixelAReverse = real(source(128).pixelAReverse);
pixelAReverse(isnan(pixelAReverse)) = 0;
pixelAReverse = pixelAReverse/max(pixelAReverse(:));
surf(grid.xAxis, grid.yAxis, pixelAReverse', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction - Sensor 128');
saveas(gcf, 'Example36_RT_sensor128', 'png'); 
saveas(gcf, 'Example36_RT_sensor128.fig');

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
