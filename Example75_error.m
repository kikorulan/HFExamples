% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex75_3D_thinveins_het;

clear all;
close all;

% Functions
norm_distance = @(x, y) sqrt(sum((x(:) - y(:)).*(x(:) - y(:))));

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('./input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%========================================================================================================================
% INITIAL PRESSURE
%========================================================================================================================
% Load Initial Pressure
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
a = axes;
t = title('Initial Pressure');
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, './figures/Example75_initial_pressure.fig');
%saveas(gcf, './figures/Example75_initial_pressure', 'epsc');

% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_14400sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);
figure();
imagesc(y0);
box on;
colorbar();
%saveas(gcf, './figures/Example75_forward_signal.fig');
%saveas(gcf, './figures/Example75_forward_signal', 'epsc');

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 5;
%============================================================
% PARAMETERS
%============================================================
% GD *************************************
GD = [];
GD.tau    = '8e18';
GD.lambda = '1e-2';
GD.iter   = int2str(iter);
% S-GD ***********************************
SGD = [];
SGD.tau    = '8e18';
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
PDHG.sigma  = '1';
PDHG.tau    = '2e18';
PDHG.theta  = '1';
PDHG.lambda = '1e-2';
PDHG.iter   = int2str(iter);
% S-PDHG *********************************
SPDHG = [];
SPDHG.sigma  = '5e-1';
SPDHG.tau    = '1.5e19';
SPDHG.theta  = '1';
SPDHG.lambda = '3e-4';
SPDHG.batch  = '90';
SPDHG.epoch  = int2str(iter);

%==============================
% Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - heterogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/Example75_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);
saveas(gcf, ['./figures/Example75_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter], 'epsc');

%==============================
% Stochastic Gradient Descent
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - heterogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/Example75_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);
saveas(gcf, ['./figures/Example75_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch], 'epsc');

%==============================
% FISTA
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FISTA - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - heterogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/Example75_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);
saveas(gcf, ['./figures/Example75_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter], 'epsc');

%==============================
% PDHG
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - heterogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/Example75_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);
saveas(gcf, ['./figures/Example75_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter], 'epsc');

%==============================
% S-PDHG
%==============================
pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - heterogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, ['./figures/Example75_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);
saveas(gcf, ['./figures/Example75_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch], 'epsc');


%========================================================================================================================
% PRIMAL DISTANCE ERROR
%========================================================================================================================
disp('******* PRIMAL DISTANCE ********');
% Gradient descent
disp('GD');
GD.tau = '2e18';
GD_error_pd = norm_distance(u0, 0*u0);
for iter = 1:5
    ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    GD_error_pd = [GD_error_pd norm_distance(u0, pp)];
end

% Stochastic Gradient descent
disp('S-GD');
SGD_error_pd = norm_distance(u0, 0*u0);
for iter = 1:5
    ppmatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    SGD_error_pd = [SGD_error_pd norm_distance(u0, pp)];
end
%%  
%%  % FISTA
%%  disp('FISTA');
%%  FISTA_error_pd = norm_distance(u0, 0*u0);
%%  for iter = 1:5
%%      ppmatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%      pp = max(0, matrix2cube(ppmatrix, Nz));
%%      FISTA_error_pd = [FISTA_error_pd norm_distance(u0, pp)];
%%  end
%%  
%%  % PDHG
%%  disp('PDHG');
%%  PDHG_error_pd = norm_distance(u0, 0*u0);
%%  for iter = 1:5
%%      ppmatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%      pp = max(0, matrix2cube(ppmatrix, Nz));
%%      PDHG_error_pd = [PDHG_error_pd norm_distance(u0, pp)];
%%  end
%%  
%%  % SPDHG
%%  disp('S-PDHG');
%%  SPDHG_error_pd = norm_distance(u0, 0*u0);
%%  for iter = 1:5
%%      ppmatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
%%      pp = max(0, matrix2cube(ppmatrix, Nz));
%%      SPDHG_error_pd = [SPDHG_error_pd norm_distance(u0, pp)];
%%  end

% Plot
x_axis = [0 1 2 3 4 5];
figure();
semilogy(x_axis, GD_error_pd, 'Color', 'r', 'Linewidth', 1.5);
hold on;
semilogy(x_axis, SGD_error_pd, 'Color', 'g', 'Linewidth', 1.5);
%semilogy(x_axis, FISTA_error_pd, 'Color', 'b', 'Linewidth', 1.5);
%semilogy(x_axis, PDHG_error_pd, 'Color', 'm', 'Linewidth', 1.5);
%semilogy(x_axis, SPDHG_error_pd, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
%legend('GD', 'S-GD', 'FISTA', 'S-PDHG');
title('Primal Distance Error - heterogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%saveas(gcf, ['./figures/Example75_pd_error.fig']);
%saveas(gcf, ['./figures/Example75_pd_error'], 'epsc');


%========================================================================================================================
% DUAL DISTANCE ERROR
%========================================================================================================================
disp('******* DUAL DISTANCE ********');
% Gradient descent
disp('GD');
GD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    GD_error_dd = [GD_error_dd norm_distance(y0, yi)];
end

% Stochastic Gradient descent
disp('S-GD');
SGD_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    SGD_error_dd = [SGD_error_dd norm_distance(y0, yi)];
end

% FISTA
disp('FISTA');
FISTA_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    FISTA_error_dd = [FISTA_error_dd norm_distance(y0, yi)];
end

% PDHG
disp('PDHG');
PDHG_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    yi = tSignal(2:end, :);
    PDHG_error_dd = [PDHG_error_dd norm_distance(y0, yi)];
end

% SPDHG
disp('S-PDHG');
SPDHG_error_dd = norm_distance(y0, 0*y0);
for iter = 1:5
    tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', int2str(iter), '.dat'], ' ', 0); 
    yi = tSignal(2:end, :);
    SPDHG_error_dd = [SPDHG_error_dd norm_distance(y0, yi)];
end

% Plot
x_axis = [0 1 2 3 4 5];
figure();
semilogy(x_axis, GD_error_dd, 'Color', 'r', 'Linewidth', 1.5);
hold on;
semilogy(x_axis, SGD_error_dd, 'Color', 'g', 'Linewidth', 1.5);
semilogy(x_axis, FISTA_error_dd, 'Color', 'b', 'Linewidth', 1.5);
semilogy(x_axis, PDHG_error_dd, 'Color', 'm', 'Linewidth', 1.5);
semilogy(x_axis, SPDHG_error_dd, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Dual Distance Error - heterogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
saveas(gcf, ['./figures/Example75_dd_error.fig']);
saveas(gcf, ['./figures/Example75_dd_error'], 'epsc');
