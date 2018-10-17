

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;
%cd /home/kiko/Documents/HighFreqCode/Examples/Ex51_reconstruction2D;
close all;
clear all;

drawForward = 1;
drawAdjoint = 1;
drawMix = 1;


%===============================================================================================================
%===============================================================================================================
%=====================                 FORWARD PROBLEM              ============================================
%===============================================================================================================
%===============================================================================================================
if(drawForward)

load gridRT_impulse;
load signalRT_nonsmooth;
load signalRT_revA_nonsmooth;
load sensor_data;


cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
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

% RT revA nonsmooth
unsorted_data = signalRT_revA_nonsmooth;
sort_data = unsorted_data;
for i = 1:Ny-2
    sort_data(Nx + i, :) = unsorted_data(Nx + 2*i, :);
    sort_data(Nx + Ny - 1 + Nx - 1 + i, :) = unsorted_data(Nx + 2*(Ny - 1) - 2*i + 1, :);
end
for j = 1:Nx
    sort_data(Nx + Ny - 2 + j, :) = unsorted_data(Nx + 2*Ny + Nx - 3 - j, :);
end
signalRT_revA_nonsmooth = sort_data; 

%========================================
% NORMS
%========================================
% k-Wave
maxKW = max(sensor_data(:));
normKW = maxKW; 
signalKW_norm = sensor_data/normKW;
% RT nonsmooth
maxRT_nonsmooth = max(signalRT_nonsmooth(:));
normRT_nonsmooth = maxRT_nonsmooth;
signalRT_nonsmooth_norm = signalRT_nonsmooth/normRT_nonsmooth;
% RT nonsmooth
maxRT_nonsmooth = max(signalRT_revA_nonsmooth(:));
normRT_nonsmooth = maxRT_nonsmooth;
signalRT_revA_nonsmooth_norm = signalRT_revA_nonsmooth/normRT_nonsmooth;


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


% RT nonsmooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalRT_nonsmooth_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.4 1]);
xlabel('t [\mus]');
%ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYNoBar);
saveas(gcf, 'Example51_RT_f2D_nonsmooth.fig');
saveas(gcf, 'Example51_RT_f2D_nonsmooth', 'png');


% RT nonsmooth
figure;
surf(1e6*Rgrid.tForward, 1:nSources, signalRT_revA_nonsmooth_norm, 'EdgeColor', 'none');
axis([0 40 1 nSources]);
view(2);
box on;
caxis([-0.4 1]);
xlabel('t [\mus]');
%ylabel('Sensor');
set(gca,'FontSize',fontSize);
set(gcf, 'pos', positionNoYNoBar);
saveas(gcf, 'Example51_RT_revA_f2D_smooth.fig');
saveas(gcf, 'Example51_RT_revA_f2D_smooth', 'png');

end
