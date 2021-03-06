

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;
%cd /home/kiko/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;
close all;

drawForward = 0;
drawAdjoint = 1;
drawMix = 1;


%===============================================================================================================
%===============================================================================================================
%=====================                 FORWARD PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawForward)

load gridRT_impulse;
load signalRT_smooth;
load signalRT_nonsmooth;
load signalGB_smooth;
load signalGB_nonsmooth;
load sensor_data;

%===============================================================================================================
% SINOGRAM
%===============================================================================================================
nSources = 764;
dcol_pos = @(x) [x(:, end) x(:, 1:end-1)];
dcol_neg = @(x) [x(:, 2:end) x(:, 1)];

positionYNoBar     = [700 700 550 630];
positionNoYBar     = [700 700 600 630];
positionYBar       = [700 700 620 630];
positionNoYNoBar   = [700 700 530 630];
%positionNoYBar = [700 700 610 630];
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% SORT DATA
%========================================
Nx = 128;
Ny = 256;
% k-Wave 
unsorted_data = sensor_data;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
sensor_data = sort_data; 

% RT smooth
unsorted_data = signalRT_smooth;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
signalRT_smooth = sort_data; 

% RT nonsmooth
unsorted_data = signalRT_nonsmooth;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
signalRT_nonsmooth = sort_data; 

% GB smooth
unsorted_data = signalGB_smooth;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
signalGB_smooth = sort_data; 

% GB nonsmooth
unsorted_data = signalGB_nonsmooth;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
signalGB_nonsmooth = sort_data; 

%========================================
% NORMS
%========================================
% k-Wave
maxKW = max(sensor_data(:));
normKW = maxKW; 
signalKW_norm = sensor_data/normKW;
% RT smooth
maxRT_smooth = max(signalRT_smooth(:));
normRT_smooth = maxRT_smooth;
signalRT_smooth_norm = signalRT_smooth/normRT_smooth;
% RT nonsmooth
maxRT_nonsmooth = max(signalRT_nonsmooth(:));
normRT_nonsmooth = maxRT_nonsmooth;
signalRT_nonsmooth_norm = signalRT_nonsmooth/normRT_nonsmooth;
% GB smooth
maxGB_smooth = max(signalGB_smooth(:));
normGB_smooth = maxGB_smooth;
signalGB_smooth_norm = real(signalGB_smooth/normGB_smooth);
% GB nonsmooth
maxGB_nonsmooth = max(signalGB_nonsmooth(:));
normGB_nonsmooth = maxGB_nonsmooth;
signalGB_nonsmooth_norm = real(signalGB_nonsmooth/normGB_nonsmooth);

%========================================
% FIGURES - FORWARD
%========================================
fontSize = 16;
% k-Wave
figure;
surf(1e6*kgrid.t_array, 1:nSources, signalKW_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
%colorbar();
caxis([-0.4 1]);
xlabel('t [\mus]');
ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionYNoBar);
saveas(gcf, 'Example51_kWave_f2D.fig');
saveas(gcf, 'Example51_kWave_f2D', 'png');


% RT smooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalRT_smooth_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.4 1]);
xlabel('t [\mus]');
%ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYNoBar);
saveas(gcf, 'Example51_RT_f2D_smooth.fig');
saveas(gcf, 'Example51_RT_f2D_smooth', 'png');
% RT nonsmooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalRT_nonsmooth_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
colorbar();
caxis([-0.4 1]);
xlabel('t [\mus]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYBar);
saveas(gcf, 'Example51_RT_f2D_nonsmooth.fig');
saveas(gcf, 'Example51_RT_f2D_nonsmooth', 'png');
% GB smooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalGB_smooth_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.4 1]);
xlabel('t [\mus]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYNoBar);
saveas(gcf, 'Example51_GB_f2D_smooth.fig');
saveas(gcf, 'Example51_GB_f2D_smooth', 'png');
% GB nonsmooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalGB_nonsmooth_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
colorbar();
caxis([-0.4 1]);
xlabel('t [\mus]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYBar);
saveas(gcf, 'Example51_GB_f2D_nonsmooth.fig');
saveas(gcf, 'Example51_GB_f2D_nonsmooth', 'png');

%========================================
% FIGURES - ERROR
%========================================
% Error
signal = dcol_neg(signalRT_smooth_norm);
error_RTsmooth = signal - signalKW_norm;

% Error - RT smooth
signal = dcol_neg(signalRT_smooth_norm);
error_RTsmooth = signal - signalKW_norm;
figure;
surf(1e6*Rgrid.tForward, 1:nSources, error_RTsmooth, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.5 0.5]);
xlabel('t [\mus]');
ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionYNoBar);
saveas(gcf, 'Example51_error_RT_f2D_smooth.fig');
saveas(gcf, 'Example51_error_RT_f2D_smooth', 'png');
% RT nonsmooth
signal = dcol_neg(signalRT_nonsmooth_norm);
error_RTnonsmooth = signal - signalKW_norm;
figure;
surf(1e6*Rgrid.tForward, 1:nSources, error_RTnonsmooth, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
colorbar();
caxis([-0.5 0.5]);
xlabel('t [\mus]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYBar);
saveas(gcf, 'Example51_error_RT_f2D_nonsmooth.fig');
saveas(gcf, 'Example51_error_RT_f2D_nonsmooth', 'png');
% GB smooth
signal = dcol_neg(signalGB_smooth_norm);
error_GBsmooth = signal - signalKW_norm;
figure;
surf(1e6*Rgrid.tForward, 1:nSources, error_GBsmooth, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.5 0.5]);
xlabel('t [\mus]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYNoBar);
saveas(gcf, 'Example51_error_GB_f2D_smooth.fig');
saveas(gcf, 'Example51_error_GB_f2D_smooth', 'png');
% GB nonsmooth
signal = dcol_neg(signalGB_nonsmooth_norm);
error_GBnonsmooth = signal - signalKW_norm;
figure;
surf(1e6*Rgrid.tForward, 1:nSources, error_GBnonsmooth, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.5 0.5]);
colorbar();
xlabel('t [\mus]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYBar);
saveas(gcf, 'Example51_error_GB_f2D_nonsmooth.fig');
saveas(gcf, 'Example51_error_GB_f2D_nonsmooth', 'png');

end

%===============================================================================================================
%===============================================================================================================
%=====================                 ADJOINT PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawAdjoint)

load recon_data_adjoint.mat;
load gridRT_impulse;
load adjointRT_smooth;
load adjointRT_nonsmooth;
load adjointGB_smooth;
load adjointGB_nonsmooth;

%========================================
% Draw parameters
%========================================
axisGrid = [0 1e3*(Rgrid.Nx-1)*Rgrid.dx 0 1e3*(Rgrid.Ny-1)*Rgrid.dy];
position     = [700 700 320 630];
positionY    = [700 700 340 630];
positionBar  = [700 700 380 630];
positionYBar = [700 700 410 630];
fontSize = 15;
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% kWave Adjoint
%========================================
pixelAReverse = p0_recon_adjoint.p_final;
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
saveas(gcf, 'Example51_kWave_recon', 'png');
saveas(gcf, 'Example51_kWave_recon.fig');

%========================================
% kWave Adjoint
%========================================
pixelAReverse = p0_recon_adjoint.p_final;
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example51_kWave_recon_2', 'png');
saveas(gcf, 'Example51_kWave_recon_2.fig');

%========================================
% RT nonsmooth
%========================================
pixelAReverse = adjointRT_nonsmooth;
maxPixelRT = max(real(pixelAReverse(:)));
pixelRT = max(0, real(pixelAReverse)/maxPixelRT);
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example51_RT_recon_nonsmooth', 'png');
saveas(gcf, 'Example51_RT_recon_nonsmooth.fig');

%========================================
% GB nonsmooth
%========================================
pixelAReverse = adjointGB_nonsmooth;
maxPixelRT = max(real(pixelAReverse(:)));
pixelRT = max(0, real(pixelAReverse)/maxPixelRT);
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example51_GB_recon_nonsmooth', 'png');
saveas(gcf, 'Example51_GB_recon_nonsmooth.fig');

%========================================
% Error RT
%========================================
adjointKW = max(0, p0_recon_adjoint.p_final/max(p0_recon_adjoint.p_final(:)));
adjointRT = max(0, adjointRT_nonsmooth/max(adjointRT_nonsmooth(:)));
pixelKWave = adjointKW - adjointRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
caxis([-.5 .5]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example51_error_RT', 'png');
saveas(gcf, 'Example51_error_RT.fig');

%========================================
% Error GB
%========================================
adjointKW = max(0, p0_recon_adjoint.p_final/max(p0_recon_adjoint.p_final(:)));
adjointRT = max(0, real(adjointGB_nonsmooth)/max(real(adjointGB_nonsmooth(:))));
pixelKWave = adjointKW - adjointRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
caxis([-.5 .5]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', position);
saveas(gcf, 'Example51_error_GB', 'png');
saveas(gcf, 'Example51_error_GB.fig');

end


%===============================================================================================================
%===============================================================================================================
%=====================             MIX ADJOINT PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawMix)
    

load recon_data_adjoint;
load adjointKWaveForward_RT;
load adjointKWaveForward_GB;
load adjointRTForward_kWave;
load adjointGBForward_kWave;
figure;
%========================================
% Draw parameters
%========================================
axisGrid = [0 1e3*(Rgrid.Nx-1)*Rgrid.dx 0 1e3*(Rgrid.Ny-1)*Rgrid.dy];
position     = [700 700 320 630];
positionY    = [700 700 340 630];
positionBar  = [700 700 380 630];
positionYBar = [700 700 410 630];
fontSize = 15;
set(0,'DefaultFigurePaperPositionMode','auto');


%========================================
% RT Adjoint - kWave data
%========================================
pixelAReverse = real(adjointKWaveForward_RT);
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example51_RT_recon_kWave_data', 'png');
saveas(gcf, 'Example51_RT_recon_kWave_data.fig');

%========================================
% GB Adjoint - kWave data
%========================================
pixelAReverse = real(adjointKWaveForward_GB);
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example51_GB_recon_kWave_data', 'png');
saveas(gcf, 'Example51_GB_recon_kWave_data.fig');

%========================================
% kWave Adjoint - RT data
%========================================
pixelAReverse = real(adjointRTForward_kWave.p_final);
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionY);
saveas(gcf, 'Example51_kWave_recon_RT_data', 'png');
saveas(gcf, 'Example51_kWave_recon_RT_data.fig');

%========================================
% kWave Adjoint - GB data
%========================================
pixelAReverse = real(adjointGBForward_kWave);
pixelKWave = max(0, pixelAReverse/max(pixelAReverse(:)));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example51_kWave_recon_GB_data', 'png');
saveas(gcf, 'Example51_kWave_recon_GB_data.fig');


%========================================
% Error mix kWave Forward - RT recon
%========================================
adjointKW = max(0, p0_recon_adjoint.p_final/max(p0_recon_adjoint.p_final(:)));
adjointRT = max(0, real(adjointKWaveForward_RT)/max(real(adjointKWaveForward_RT(:))));
pixelKWave = adjointKW - adjointRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
caxis([-.5 .5]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example51_mix_error_RT', 'png');
saveas(gcf, 'Example51_mix_error_RT.fig');

%========================================
% Error mix kWave Forward - GB recon
%========================================
adjointKW = max(0, p0_recon_adjoint.p_final/max(p0_recon_adjoint.p_final(:)));
adjointRT = max(0, real(adjointKWaveForward_GB)/max(real(adjointKWaveForward_GB(:))));
pixelKWave = adjointKW - adjointRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
caxis([-.5 .5]);
box on;
xlabel('x [mm]');
%ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionBar);
saveas(gcf, 'Example51_mix_error_GB', 'png');
saveas(gcf, 'Example51_mix_error_GB.fig');


end


