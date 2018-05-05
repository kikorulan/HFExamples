%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;

close all;
% Measure computational time
tic;
start_time = clock;

load recon_data_RT.mat;
load sensor_data.mat;
%load recon_data_sensors.mat;
load recon_data_adjoint.mat;

%==================================================================================
% SIMULATION
%==================================================================================
%==============================
% RT
%==============================
% Define source
source_sub_128 = source(1:2:end);
source_sub_64 = source(2:4:end);
source_sub_32 = source(4:8:end);
source_sub_16 = source(8:16:end);

% Choose source
sourceRT = source_sub_16;

% Define parameters for 16 sensors
lambda = 3e-12;
lipschitz = 1.5e-10;
nIter = 5;

%%  % Define parameters for 64 sensors
%%  lambda = 1e-11;
%%  lipschitz = 5e-10;
%%  nIter = 5;

% Define operators
norm_factor = max(grid.u(:))/max(pixelAReverse_16(:));
%norm_factor = 1;
forward_operator = @(initial_pressure) grid.operator_forward(sourceRT, initial_pressure, 1);
inverse_operator = @(forward_data) grid.operator_inverse(sourceRT, forward_data);

% Initialise forward and inverse data
%clear b;
for n = 1:length(sourceRT)
    if(isempty(sourceRT(n).aForward_initial))
        sourceRT(n).setForwardSignal_initial(sourceRT(n).aForward);
    end
    b(n, :) = sourceRT(n).aForward_initial;
end
x0 = pixelAReverse_16;

% Run iterations
%pixelRT = fista(b, norm_factor*x0, lambda, lipschitz, nIter, forward_operator, inverse_operator);
%pixelRT = gradient_descent(b, x0, 1, nIter, forward_operator, inverse_operator);

%==============================
% k-Wave
%==============================
% Define subsampled data
sensor_data_128 = sensor_data(1:2:end, :);
sensor_data_64 = sensor_data(2:4:end, :);
sensor_data_32 = sensor_data(4:8:end, :);
sensor_data_16 = sensor_data(8:16:end, :);

% Define the sensors
sensor_128.mask = zeros(grid.Nx, grid.Ny);
sensor_128.mask(1, 1:2:end) = 1;
sensor_64.mask = zeros(grid.Nx, grid.Ny);
sensor_64.mask(1, 2:4:end) = 1;
sensor_32.mask = zeros(grid.Nx, grid.Ny);
sensor_32.mask(1, 4:8:end) = 1;
sensor_16.mask = zeros(grid.Nx, grid.Ny);
sensor_16.mask(1, 8:16:end) = 1;

% Choose sensor
sensor_dataKW = sensor_data_16;
sensorKW = sensor_16;
p0_reconKW = p0_recon_adjoint_16;

% Define parameters for 16 sensors
lambda = 8e-3;
lipschitz = 4e-1;
nIter = 5;

%%  % Define parameters with 64 sensors
%%  lambda = 4e-2;
%%  lipschitz = 2;
%%  nIter = 5;

% Define operators
forward_operator = @(initial_pressure) forward_operator_kWave(kgrid, medium, initial_pressure, sensorKW, input_args);
inverse_operator = @(forward_data) inverse_operator_kWave(kgrid, medium, forward_data, sensorKW, input_args);

% Initialise forward and inverse data
b = sensor_dataKW;
x0 = p0_reconKW.p_final;
pixelKWave_ini = x0;

% Run iterations
pixelKWave = fista(b, x0, lambda, lipschitz, nIter, forward_operator, inverse_operator);
%==================================================================================
% VISUALISATION
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;
axisGrid = [0 grid.xAxis(end) 0 grid.yAxis(end)];

position = [700 700 320 600];
positionNoY = [700 700 300 600];
positionNoYBar = [700 700 363 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Initial Pressure
%==============================
figure;
surf(grid.xAxis, grid.yAxis, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Initial Pressure');

%==============================
% Reconstruction - RT initial
%==============================
pixelRT_ini = max(0, pixelAReverse_16)/max(pixelAReverse_16(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelRT_ini', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
%set(gca, 'YTick', []);
set(gcf, 'pos', position);
%title('Reconstruction RT');
saveas(gcf, 'Example33_RT_16', 'png');

%==============================
% Reconstruction - RT FISTA
%==============================
pixelRT = max(0, pixelRT)/max(pixelRT(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
%set(gca, 'YTick', []);
set(gcf, 'pos', position);
%title('Reconstruction RT');
saveas(gcf, 'Example33_RT_16_FISTA', 'png');

%========================================
% Reconstruction - kWave
%========================================
pixelKWave_ini = max(0, pixelKWave_ini)/max(pixelKWave_ini(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelKWave_ini', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
%title('Reconstruction k-Wave (TR)');
saveas(gcf, 'Example33_kWave_16', 'png');
saveas(gcf, 'Example33_kWave_16.fig');

%========================================
% Reconstruction - kWave FISTA
%========================================
pixelKWave = max(0, pixelKWave)/max(pixelKWave(:));
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
%title('Reconstruction k-Wave (TR)');
saveas(gcf, 'Example33_kWave_16_FISTA', 'png');
saveas(gcf, 'Example33_kWave_16_FISTA.fig');

%========================================
% Error - RT initial
%========================================
errorRT_TR = pixelRT_ini - pixelKWave_ini;
%errorRT_TR = grid.u - pixelKWave_ini;
figure;
surf(grid.xAxis, grid.yAxis, errorRT_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
saveas(gcf, 'Example33_RT_error_16', 'png'); 

%========================================
% Error - RT FISTA
%========================================
errorRT_TR = pixelRT - pixelKWave;
%errorRT_TR = grid.u - pixelKWave;
figure;
surf(grid.xAxis, grid.yAxis, errorRT_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
%title('Error RT FISTA');
saveas(gcf, 'Example33_RT_error_16_FISTA', 'png'); 

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

