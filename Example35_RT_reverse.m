%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;

close all;
%clear all;

%load gridRT.mat;
load recon_data.mat;
% Measure computational time
tic;
start_time = clock;
 
%run colourMap;
 
%========================================
% Compute filters
%========================================
grid.timeReverseFilter(100);

%========================================
% Compute reverse signal
%========================================

for n = 1:(256+127)
    grid.timeReverseBeam(n);
end
%grid.computeTimeReverseParallel();

pixelAReverseSensors = grid.timeReverseSignal();


%save gridRT_recon.mat grid -v7.3;
%================================================================================
% VISUALISATION
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex35_reconstruction;

X = 0:grid.dx:(grid.Nx-1)*grid.dx;
Y = 0:grid.dy:(grid.Ny-1)*grid.dy;
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
saveas(gcf, 'Example35_C', 'png'); 
saveas(gcf, 'Example35_C.fig'); 

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
surf(X, Y, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Initial Pressure');
saveas(gcf, 'Example35_U', 'png'); 
saveas(gcf, 'Example35_U.fig'); 

%========================================
% Reconstruction - RT
%========================================
pixelRT = real(grid.pixelAReverse);
figure;
surf(X, Y, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction - RT');
saveas(gcf, 'Example35_RT_recon', 'png');
saveas(gcf, 'Example35_RT_recon.fig');

%========================================
% Reconstruction - kWave
%========================================
pixelKWave = max(0, p0_recon)/max(p0_recon(:));
figure;
surf(X, Y, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction - kWave');
saveas(gcf, 'Example35_kWave_recon', 'png');
saveas(gcf, 'Example35_kWave_recon.fig');

%========================================
% Error - RT-kWave
%========================================
errorRT = pixelKWave - pixelRT;
figure;
surf(X, Y, errorRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Error: RT-kWave');
saveas(gcf, 'Example35_RT_error', 'png');
saveas(gcf, 'Example35_RT_error.fig');

%========================================
% Error - kWave
%========================================
errorKWave = grid.u - pixelKWave;
figure;
surf(X, Y, errorKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-1, 1]);
colorbar();
title('Error: kWave');
saveas(gcf, 'Example35_kWave_error.fig');

%========================================
% Reconstruction - Sensor 1
%========================================
figure;
pixelAReverse = grid.source(128).pixelAReverse;
pixelAReverse(isnan(pixelAReverse)) = 0;
pixelAReverse = pixelAReverse/max(pixelAReverse(:));
surf(X, Y, pixelAReverse', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction - Sensor 128');
saveas(gcf, 'Example35_RT_sensor128', 'png'); 
saveas(gcf, 'Example35_RT_sensor128.fig');

%%  %========================================
%%  % Angle Correction - Sensor 1
%%  %========================================
%%  figure;
%%  pixelAngle = reshape(grid.source(128).pixelAngleCorrection, 128, 256);
%%  %pixelAngle(isnan(pixelAngle)) = 0;
%%  surf(X, Y, pixelAngle', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  colorbar();
%%  box on;
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  set(gcf, 'pos', position);
%%  title('Angle Correction - Sensor 128');
%%  saveas(gcf, 'Example35_RT_sensor128_angle', 'png'); 
%%  saveas(gcf, 'Example35_RT_sensor128_angle.fig');

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
