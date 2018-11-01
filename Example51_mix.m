%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex51_reconstruction2D;

close all;
%clear all;

%load RgridRT.mat;
load signalRT_nonsmooth.mat
load recon_data_adjoint.mat
%load recon_data.mat;
load sensor_data.mat
% Measure computational time
tic;
start_time = clock;
 
%run colourMap;

%================================================================================
% kWave reconstruction using RT data
%================================================================================

sensor_data_RT = signalRT_nonsmooth;
% run the time reversal reconstruction
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% Build sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny);
sensor_adjoint.record = {'p_final'};
source_adjoint.p_mask = zeros(kgrid.Nx, kgrid.Ny);
source_adjoint.p_mask(1, :) = 1;
source_adjoint.p_mask(end, :) = 1;
source_adjoint.p_mask(:, 1) = 1;
source_adjoint.p_mask(:, end) = 1;
source_adjoint.p = fliplr(sensor_data);
p0_recon_adjoint = kspaceFirstOrder2D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:});

%================================================================================
% RT reconstruction using kWave data
%================================================================================
%%  %========================================
%%  % Compute filters
%%  %========================================
%%  Rgrid.inverse_filter(100);
%%  
%%  %========================================
%%  % Compute reverse signal
%%  %========================================
%%  nSources = 764;
%%  for n = 1:nSources
%%      source(n).setForwardSignal(sensor_data(n, :));
%%  end
%%  for n = 1:nSources
%%      disp(n)
%%      Rgrid.inverse_beam(source(n));
%%  end
%%  %Rgrid.computeAdjointParallel(source);
%%  
%%  pixelAReverseSensors = Rgrid.inverse_signal(source);


% adjointKWaveForward_RT = Rgrid.pixelAReverse;
adjointRTForward_kWave = p0_recon_adjoint;
% save adjointKWaveForward_RT adjointKWaveForward_RT;
save adjointRTForward_kWave adjointRTForward_kWave;

%save RgridRT_recon.mat Rgrid -v7.3;
%================================================================================
% VISUALISATION
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;

axisGrid = [0 1e3*Rgrid.xAxis(end) 0 1e3*Rgrid.yAxis(end)];

position     = [700 700 300 630];
positionY    = [700 700 320 630];
positionBar  = [700 700 363 630];
positionYBar = [700 700 390 630];
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% Reconstruction - RT
%========================================
%%  % Using k-Wave sensors
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
%%  saveas(gcf, 'Example51_mix_nonsmoothGB-kWave_recon', 'png');
%%  saveas(gcf, 'Example51_mix_nonsmoothGB-kWave_recon.fig');

%========================================
% Reconstruction - kWave Adjoint
%========================================
%%  % Using RT sensors
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
%%  saveas(gcf, 'Example51_mix_kWave-GB_recon', 'png');
%%  saveas(gcf, 'Example51_mix_kWave-GB_recon.fig');

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
%%  saveas(gcf, 'Example51_mix_error_GB_nonsmooth','png');
%%  saveas(gcf, 'Example51_mix_error_GB_nonsmooth.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

