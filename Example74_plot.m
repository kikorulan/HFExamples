% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex74_3D_thinveins;

clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;


%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%============================================================================================================================================
% FORWARD PROBLEM
%============================================================================================================================================

%==================================================
% TIME SIGNAL - RT
%==================================================
% Import data
filenameData = './output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

% Import data
filenameData = './input_data/forwardSignal_reference_14400sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT2 = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT2);
box on;
colorbar();


% Difference
figure
imagesc(inputRT-inputRT2);
box on;
colorbar();

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('output_data/Example74_forward_output_1600sensors.h5', '/p');
inputKWave = sensor_data(:, :);
% Plot
figure;
imagesc(inputKWave);
box on;
colorbar();

%==================================================
% COMPARISON IMAGE
%==================================================
dcol_pos = @(x, n) [x(:, end-(n-1):end) x(:, 1:end-n)];
dcol_neg = @(x, n) [x(:, n+1:end) x(:, 1:n)];
%timeKWave = kgrid.t_array;
nSensors = size(sensor_data, 1);
normRT = max(inputRT(1,:));
normKWave = max(inputKWave(1,:));

for n = 1:10
    % Plot
    figure;
    imagesc(inputKWave/normKWave - dcol_pos(inputRT, n)/normRT);
    box on;
    colorbar();
end

for n = 1:10
    figure;
    imagesc(inputKWave/normKWave - dcol_neg(inputRT, n)/normRT);
    box on;
    colorbar();
end

%==================================================
% COMPARISON WITH K-WAVE
%==================================================
% Import data
sensor_data = h5read('output_data/Example74_forward_output_1600sensors.h5', '/p');
filenameData = 'output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, ' ', 0);
% Input
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
%timeKWave = kgrid.t_array;
inputKWave = sensor_data;
nSensors = size(sensor_data, 1);

normRT = max(inputRT(1,:));
normKWave = max(inputKWave(1,:));

% Norm difference
difNorm = norm(inputRT(:)/normRT - inputKWave(:)/normKWave);

% Plot all sensors
vecS = 1:100:1600;
for i = 1:length(vecS)
    % Normalisation - RT
    signalRT = inputRT(vecS(i), :)/normRT;
    % Normalisation - kWave
    signalKWave = inputKWave(vecS(i), :)/normKWave;
    % Plot
    figure;
    plot(timeRT, signalRT, 'Color', 'r', 'LineWidth', 2);
    hold on;
    set(gca,'FontSize',18);
    plot(timeRT, signalKWave, 'Color', 'blue', 'LineWidth', 2);
    axis([0 1.5e-5 -1.8 1.8]);
    legend('RT', 'k-Wave');
    xlabel('time (s)');
    ylabel('amplitude');
    ax = gca;
    ax.GridAlpha = 0.5;
    grid on;
    title(['RT vs kWave - ', int2str(i)]);
end

for i = 1:length(vecS)
    % Error
    signalRT = inputRT(vecS(i), :)/normRT;
    signalKWave = inputKWave(vecS(i), :)/normKWave;
    % Spline
    signalKWave_spline = spline(timeRT, signalKWave, timeRT);
    error = signalRT - signalKWave_spline;
    figure;
    plot(timeRT, error, 'Color', 'r', 'LineWidth', 2);
    set(gca,'FontSize',18);
    axis([0 1.5e-5 -.2 .2]);
    legend('error');
    xlabel('time (s)');
    ylabel('error');
    ax = gca;
    ax.GridAlpha = 0.5;
    grid on;
    title(['Error RT vs kWave - ', int2str(i)]);
end

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================
% Position
position     = [700 700 300 630];
positionY    = [700 700 320 630];
positionBar  = [700 700 363 630];
positionYBar = [700 700 390 630];
set(0,'DefaultFigurePaperPositionMode','auto');

%==================================================
% INITIAL PRESSURE
%==================================================
% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
a = axes;
t = title('Initial Pressure');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example74_initial_pressure.fig');
% Load Initial Pressure
load ./input_data/initial_pressure_veins_smooth;
plot_projection(initial_pressure_veins_smooth, dx);

%==================================================
% Reconstruction - kWave
%==================================================
PML_size = 10;
% 14400 sensors
p0_recon_PML = h5read('output_data/Example74_adjoint_output_14400sensors.h5', '/p_final');
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
p0_recon = p0_recon/max(p0_recon(:));
plot_projection(p0_recon, dx);
a = axes;
t = title('k-Wave adjoint projection');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example74_kWave_adjoint_14400sensors.fig');

% 1600 sensors
p0_recon_PML = h5read('output_data/Example74_adjoint_output_1600sensors.h5', '/p_final');
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
p0_recon = p0_recon/max(p0_recon(:));
plot_projection(p0_recon, dx);
a = axes;
t = title('k-Wave adjoint projection');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example74_kWave_adjoint_1600sensors.fig');


%==================================================
% ADJOINT k-Wave
%==================================================
pixelPressureMatrix = importdata('input_data/pressure_adjoint_kWave_14400sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_kWave_1600sensors.dat', ' ', 0);
pixelPressure_adjoint = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure_adjoint, dx);

%==================================================
% ADJOINT RT
%==================================================
% 1600
pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_RT_1600sensors.dat', ' ', 0);
pixelPressure = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure, dx);

% 14400
pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_RT_14400sensors.dat', ' ', 0);
normRT = 5e17;
pixelPressureMatrix_norm = pixelPressureMatrix*normRT;
dlmwrite('input_data/pixelPressure_adjoint_RT_14400sensors_norm.dat', pixelPressureMatrix_norm, 'delimiter', ' ');
pixelPressure = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure, dx);
pixelPressureMatrix_norm = pixelPressureMatrix*0;
dlmwrite('input_data/pixelPressure_0.dat', pixelPressureMatrix_norm, 'delimiter', ' ');

% 14400 norm
pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_RT_14400sensors_norm.dat', ' ', 0);
pixelPressure = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure, dx);

% 0
pixelPressureMatrix = importdata('input_data/pixelPressure_0.dat', ' ', 0);
pixelPressure = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure, dx);


%==================================================
% AUXILIAR RECONSTRUCTION
%==================================================
% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_PDHG_sigma2e1_tau4e17_theta1_lambda1e-2_iter0.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_S-GD_14.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_FISTA_4.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_S-FISTA_29.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_PDHG_1.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_S-GD_tau4e18_lambda3e-4_batch90_epoch5.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_TVdenoised_1e-14-3000.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);


%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 1;
%============================================================
% PARAMETERS
%============================================================
% GD *************************************
GD = [];
GD.tau    = '1e18';
GD.lambda = '1e-2';
GD.iter   = int2str(iter);
% S-GD ***********************************
SGD = [];
SGD.tau    = '4e18';
SGD.lambda = '3e-4';
SGD.batch  = '90';
SGD.epoch  = int2str(iter);
% FISTA **********************************
FISTA = [];
FISTA.tau    = '1e18';
FISTA.lambda = '1e-2';
FISTA.iter   = int2str(iter);
% PDHG ***********************************
PDHG = [];
PDHG.sigma  = '1e0';
PDHG.tau    = '1e18';
PDHG.theta  = '1';
PDHG.lambda = '1e-2';
PDHG.iter   = int2str(iter);
% S-PDHG *********************************
SPDHG = [];
SPDHG.sigma  = '1e1';
SPDHG.tau    = '4e18';
SPDHG.theta  = '1';
SPDHG.lambda = '3e-4';
SPDHG.batch  = '90';
SPDHG.epoch  = int2str(iter);

%==============================
% Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example74_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);

%==============================
% Stochastic Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example74_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FISTA - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example74_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);

%==============================
% PDHG
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example74_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);

%==============================
% S-PDHG
%==============================
pixelPressureMatrix = importdata(['./results/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example74_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);


%========================================================================================================================
% PRIMAL DISTANCE ERROR
%========================================================================================================================

