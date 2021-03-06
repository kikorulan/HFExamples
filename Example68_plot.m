% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex68_3D_veins_resize;

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
filenameData = 'output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

% Import data
filenameData = 'input_data/forwardSignal_reference_1600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
box on;
colorbar();

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('output_data/Example68_forward_output_1600sensors.h5', '/p');
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
%load input_data/sensor_data_4balls;
sensor_data = h5read('output_data/Example68_forward_output_1600sensors.h5', '/p');
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


%==================================================
% TIME SIGNAL - RT
%==================================================
nS = 1;
% Import data - 1
timeSignal_dt1 = importdata('input_data/forwardSignal_1600sensors_dt1.6e-8_lowres_interp.dat', ' ', 0);
timeRT_dt1 = timeSignal_dt1(1, :);
inputRT_dt1 = timeSignal_dt1(2:end, :);
normDT1 = max(inputRT_dt1(nS, :));
% Import data - 2
timeSignal_dt2 = importdata('input_data/forwardSignal_1600sensors_dt1.6e-8_lowres_interp_delay2.dat', ' ', 0);
timeRT_dt2 = timeSignal_dt2(1, :);
inputRT_dt2 = timeSignal_dt2(2:end, :);
normDT2 = max(inputRT_dt2(nS, :));
%%  % Import data - 3
%%  timeSignal_dt3 = importdata('input_data/forwardSignal_1600sensors_dt1.6e-8_lowres.dat', ' ', 0);
%%  timeRT_dt3 = timeSignal_dt3(1, :);
%%  inputRT_dt3 = timeSignal_dt3(2:end, :);
%%  normDT3 = max(inputRT_dt3(nS, :));
%%  % k-Wave
%%  sensor_data = h5read('output_data/Example68_forward_output_1600sensors.h5', '/p');
%%  normKW = max(sensor_data(nS, :));
% Plot
figure;
hold on;
plot(timeRT_dt1, inputRT_dt1(nS, :)/normDT1, 'Color', 'r');
plot(timeRT_dt2, inputRT_dt2(nS, :)/normDT2, 'Color', 'g');
%plot(timeRT_dt3, inputRT_dt3(nS, :)/normDT3, 'Color', 'b');
%plot(timeRT_dt1, sensor_data(nS, :)/normKW, 'Color', 'm');
%legend('1.6e-8', '8e-9', '4e-9', 'kWave');
%legend('8e-9', '4e-9', 'kWave');
%legend('1.6e-8 9 sensors', '8e-9 lowres', '1.6e-8 lowres', 'kWave');

%==================================================
% Compute forward signal difference
%==================================================
% Import data
filenameData = 'input_data/forwardSignal_1600sensors_dt1.6e-8_lowres_interp_delay2.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
signalRT_1 = timeSignal(2:end, :);
% Import data
filenameData = 'input_data/forwardSignal_kWave_adjoint_1600sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
signalRT_2 = timeSignal(2:end, :);
% Create signal
differenceRT = signalRT_1 - signalRT_2;
forwardDif = [timeRT; differenceRT];
dlmwrite('input_data/forwardSignal_difference_1600sensors.dat', forwardDif, 'delimiter', ' ');

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
saveas(gcf, 'figures/Example68_initial_pressure.fig');
% Load Initial Pressure
load ./input_data/initial_pressure_veins_smooth;
%plot_projection(initial_pressure_veins_smooth, dx);

%==================================================
% Reconstruction - kWave
%==================================================
p0_recon_PML = h5read('output_data/Example68_adjoint_output_57600sensors.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
p0_recon = p0_recon/max(p0_recon(:));
plot_projection(p0_recon, dx);
a = axes;
t = title('k-Wave adjoint projection');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example68_kWave_adjoint_57600sensors.fig');
%==================================================
% Reconstruction - RT
%==================================================
% PLOT
pixelPressureMatrix = importdata('output_data/PixelPressure_delay_0_0.dat', ' ', 0);
pixelPressure_delay_0_0 = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_delay_0_0, dx);
saveas(gcf, 'figures/Example68_RT_adjoint_1600sensors_delay_0_0.fig');
% PLOT
pixelPressureMatrix = importdata('output_data/PixelPressure_delay_2_0.dat', ' ', 0);
pixelPressure_delay_2_0 = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_delay_2_0, dx);
saveas(gcf, 'figures/Example68_RT_adjoint_1600sensors_delay_2_0.fig');
% PLOT
pixelPressureMatrix = importdata('output_data/PixelPressure_delay_2_-1.dat', ' ', 0);
pixelPressure_delay_2_1 = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_delay_2_1, dx);
saveas(gcf, 'figures/Example68_RT_adjoint_1600sensors_delay_2_-1.fig');

% Line plot
posY = 70;
posZ = 70;
maxU0 = max(u0(:));
maxKW = max(p0_recon(:));
maxRT_1 = max(pixelPressure_delay_0_0(:));
maxRT_2 = max(pixelPressure_delay_2_0(:));
maxRT_3 = max(pixelPressure_delay_2_1(:));
lineU0 = u0(:, posY, posZ);
lineKW = p0_recon(:, posY, posZ);
lineRT_1 = pixelPressure_delay_0_0(:, posY, posZ);
lineRT_2 = pixelPressure_delay_2_0(:, posY, posZ);
lineRT_3 = pixelPressure_delay_2_1(:, posY, posZ);

figure;
hold on;
plot(lineU0/maxU0);
plot(lineKW/maxKW);
plot(lineRT_1/maxRT_1);
plot(lineRT_2/maxRT_2);
plot(lineRT_3/maxRT_3);
legend('u0', 'kW', 'RT delay 0, 0', 'RT delay 2, 0', 'RT delay 2, -1');
saveas(gcf, 'figures/Example68_line_xdim');

%==================================================
% ADJOINT
%==================================================
% Import data
pixelPressureMatrix = importdata('input_data/pressure_kWave_adjoint_57600sensors.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('input_data/pixelPressure_kWave_adjoint_1600sensors.dat', ' ', 0);
pixelPressure_adjoint = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(p, dx);
plot_projection(-pixelPressure_adjoint, dx);

%==================================================
% ITERATIVE RECONSTRUCTION
%==================================================
% Gradient
pixelPressureMatrix = importdata('results/pixelPressure_GD_tau1e18_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('GD - t = 1e18, l = 1e-2, iter = 50 - homogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example68_GD_tau1e18_lambda1e-2_iter50.fig');

% Stochastic gradient
pixelPressureMatrix = importdata('results/pixelPressure_S-GD_tau2e18_lambda3e-4_epoch50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('S-GD - t = 2e18, l = 3e-4, iter = 50 - homogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example68_S-GD_tau2e18_lambda3e-4_iter50.fig');

% FISTA
pixelPressureMatrix = importdata('results/pixelPressure_FISTA_tau1e18_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('FISTA - t = 1e18, l = 1e-2, iter = 50 - homogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example68_FISTA_tau1e18_lambda1e-2_iter50.fig');

% Stochastic FISTA
%pixelPressureMatrix = importdata('results/pixelPressure_S-FISTA_tau1e18_lambda5e-4_epoch50.dat', ' ', 0);
%pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%plot_projection(pixelPressure, dx);
%a = axes;
%t = title('GD - t = 1e18, l = 1e-2, iter = 50 - homogeneous SS');
%a.Visible = 'off'; 
%t.Visible = 'on'; 
%saveas(gcf, 'Example68_GD_tau1e18_lambda1e-2_iter50.fig');

% PDHG
pixelPressureMatrix = importdata('results/pixelPressure_PDHG_sigma1_tau1e18_theta1_lambda1e-2_iter50.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title('PDHG - s = 1, t = 1e18, l = 1e-2, iter = 50 - homogeneous SS');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, 'figures/Example68_PDHG_sigma1_tau1e18_lambda1e-2_iter50.fig');

%==================================================
% AUXILIAR RECONSTRUCTION
%==================================================
% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_GD_49.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_SGD_49.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_FISTA_22.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_S-FISTA_29.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_PDHG_48.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_TVdenoised_1e-14-3000.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

