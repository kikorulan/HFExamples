%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

close all;
%clear all;

%load gridRT.mat;
load recon_data.mat;
load sensor_data.mat
% Measure computational time
tic;
start_time = clock;
 
%run colourMap;

%================================================================================
% kWave reconstruction using RT data
%================================================================================

% Assign data
sensor_data_RT = aForward_RT;
% Reconstruct the initial pressure
p_xy = kspaceLineRecon(sensor_data_RT.', dy, dt, medium.sound_speed, 'Plot', true, 'PosCond', true);
% define a second k-space grid using the dimensions of p_xy
[Nx_recon, Ny_recon] = size(p_xy);
kgrid_recon = makeGrid(Nx_recon, dt*medium.sound_speed, Ny_recon, dy);
% resample p_xy to be the same size as source.p0
p_xy_rs_RT = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));

%================================================================================
% RT reconstruction using kWave data
%================================================================================
for n = 1:256
    grid.setForwardSignal(n, sensor_data(n, :));
end
%========================================
% Compute filters
%========================================
grid.timeReverseFilter(10);

%========================================
% Compute reverse signal
%========================================

for n = 1:nSources
    grid.timeReverseBeam(n);
end
pixelAReverseSensors = grid.timeReverseSignal();

%grid.computeTimeReverseParallel();
%save gridRT_recon.mat grid -v7.3;
%================================================================================
% VISUALISATION
%================================================================================
X = 0:grid.dx:(grid.Nx-1)*grid.dx;
Y = 0:grid.dy:(grid.Ny-1)*grid.dy;
axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 400 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%========================================
% Reconstruction - RT
%========================================
% Using kWave sensors
figure;
surf(X, Y, grid.pixelAReverse', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction - RT using kWave sensors');
saveas(gcf, 'Example31_mix_RT-kWave_recon','png'); 
saveas(gcf, 'Example31_mix_RT-kWave_recon.fig');

%========================================
% Reconstruction - kWave
%========================================
% Using RT sensors
pixelKWave_RT = p_xy_rs_RT/max(abs(p_xy_rs_RT(:)));
pixelKWave = p_xy_rs/max(abs(p_xy_rs(:)));
figure;
surf(X, Y, pixelKWave_RT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Reconstruction - kWave using RT sensors');
saveas(gcf, 'Example31_mix_kWave-RT_recon','png'); 
saveas(gcf, 'Example31_mix_kWave-RT_recon.fig');

%========================================
% Error - RT using kWave sensors
%========================================
errorRT = pixelKWave - grid.pixelAReverse;
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
%title('Error - RT using kWave data');
saveas(gcf, 'Example31_mix_RT-kWave_error', 'png');
saveas(gcf, 'Example31_mix_RT-kWave_error.fig');

%========================================
% Error - kWave using RT sensors
%========================================
errorKWave = pixelKWave - pixelKWave_RT;
figure;
surf(X, Y, errorKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-0.5, 0.5]);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Error - kWave using RT sensor data');
saveas(gcf, 'Example31_mix_kWave-RT_error', 'png');
saveas(gcf, 'Example31_mix_kWave-RT_error.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
