%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;

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
Rgrid.inverse_filter(100);

%========================================
% Compute reverse signal
%========================================
%Rgrid.inverse_beam(1);

%aReverse = Rgrid.inverse_beam_adjoint();
nSources = 764;
for n = 1:nSources
    disp(n)
    Rgrid.inverse_beam(source(n));
end
%Rgrid.computeAdjointParallel(source);

pixelAReverseSensors = Rgrid.inverse_signal(source);

adjointGB_nonsmooth = Rgrid.pixelAReverse;
save adjointGB_nonsmooth adjointGB_nonsmooth;

%================================================================================
% VISUALISATION
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;

axisGrid = [0 1e3*(Rgrid.Nx-1)*Rgrid.dx 0 1e3*(Rgrid.Ny-1)*Rgrid.dy];

position     = [700 700 300 630];
positionY    = [700 700 320 630];
positionBar  = [700 700 363 630];
positionYBar = [700 700 390 630];
set(0,'DefaultFigurePaperPositionMode','auto');


%========================================
% Reconstruction - RT
%========================================
%%  pixelAReverse = Rgrid.pixelAReverse;
%%  maxPixelRT = max(real(pixelAReverse(:)));
%%  pixelRT = max(0, real(pixelAReverse)/maxPixelRT);
%%  figure;
%%  surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  %colorbar();
%%  box on;
%%  xlabel('x [mm]');
%%  ylabel('y [mm]');
%%  set(gcf, 'pos', positionY);
%%  %title('Reconstruction RT');
%%  saveas(gcf, 'Example51_RT_recon_GB_nonsmooth', 'png');
%%  saveas(gcf, 'Example51_RT_recon_GB_nonsmooth.fig');

%========================================
% Reconstruction - kWave Adjoint
%========================================
%%  pixelKWave = max(0, p0_recon_adjoint.p_final/max(p0_recon_adjoint.p_final(:)));
%%  figure;
%%  surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  colorbar();
%%  box on;
%%  xlabel('x [mm]');
%%  %ylabel('y (m)');
%%  %set(gca, 'YTick', []);
%%  set(gcf, 'pos', positionYBar);
%%  %title('Reconstruction k-Wave');
%%  saveas(gcf, 'Example51_kWave_recon', 'png');
%%  saveas(gcf, 'Example51_kWave_recon.fig');

%========================================
% Error - RT-kWave
%========================================
%%  errorRT = pixelRT - pixelKWave;
%%  figure;
%%  surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, errorRT', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  caxis([-.5, .5]);
%%  colorbar();
%%  box on;
%%  xlabel('x [mm]');
%%  %ylabel('y [mm]');
%%  %set(gca, 'YTick', []);
%%  set(gcf, 'pos', positionYBar);
%%  %title('Error RT - kWave (ADJ)');
%%  saveas(gcf, 'Example51_RT_errorKW_GB_nonsmooth', 'png');
%%  saveas(gcf, 'Example51_RT_errorKW_GB_nonsmooth.fig');

%========================================
% Error - RT-Phantom
%========================================
%%  errorRT = pixelRT - Rgrid.u;
%%  figure;
%%  surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, errorRT', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  caxis([-.5, .5]);
%%  colorbar();
%%  box on;
%%  xlabel('x [mm]');
%%  %ylabel('y [mm]');
%%  set(gca, 'YTick', []);
%%  set(gcf, 'pos', positionYBar);
%%  %title('Error RT - kWave (ADJ)');
%%  saveas(gcf, 'Example51_RT_error_GB_nonsmooth', 'png');
%%  saveas(gcf, 'Example51_RT_error_GB_nonsmooth.fig');

%========================================
% Error - kWave
%========================================
%%  errorKWave = Rgrid.u - pixelKWave;
%%  figure;
%%  surf(Rgrid.xAxis, Rgrid.yAxis, errorKWave', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisRgrid);
%%  caxis([-1, 1]);
%%  colorbar();
%%  title('Error: kWave');
%%  saveas(gcf, 'Example51_kWave_error.fig');

%========================================
% Reconstruction - Sensors
%========================================
%%  meanSensors = mean(abs(pixelAReverseSensors(:)));
%%  scrollView(permute(fliplr(real(pixelAReverseSensors)), [2 1 3]), 3, [-6*meanSensors 6*meanSensors])

%========================================
% Pixel Time propagation
%========================================
%%  figure;
%%  surf(X, Y, Rgrid.pixelTime', 'EdgeColor', 'none');
%%  view(2);
%%  colorbar();
%%  title('Time Propagation');

%========================================
% Pixel Reverse Time
%========================================
%%  figure;
%%  surf(X, Y, Rgrid.pixelTReverse', 'EdgeColor', 'none');
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
