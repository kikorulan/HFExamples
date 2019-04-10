% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex81_3D_veins_subsampled_het;
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex81_3D_veins_subsampled_het;

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
% PRIMAL AND DUAL DATA
%========================================================================================================================
% Load Initial Pressure
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
a = axes;
t = title('Initial Pressure');
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, './figures/Example81_initial_pressure.fig');

% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_3600sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);

% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
h = plot_projection(pressure_adjoint, dx);

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 5;

%==============================
% Gradient Descent
%==============================
%%  % GD **************************
%%  GD = [];
%%  GD.tau    = '1.6e18';
%%  GD.lambda = '1e-3';
%%  GD.iter   = int2str(5);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['FB - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  saveas(gcf, ['./figures/Example81_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);

%==============================
% Stochastic Gradient Descent
%==============================
%%  % S-GD ************************
%%  SGD = [];
%%  SGD.tau    = '3.2e18';
%%  SGD.lambda = '1e-3';
%%  SGD.batch  = '100';
%%  SGD.epoch  = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['S-FB - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  saveas(gcf, ['./figures/Example81_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
%%  % FISTA ***********************
%%  FISTA = [];
%%  FISTA.tau    = '8e17';
%%  FISTA.lambda = '1e-3';
%%  FISTA.iter   = int2str(iter);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['AFB - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  saveas(gcf, ['./figures/Example81_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);


%==============================
% PDHG
%==============================
%%  % PDHG ************************
%%  PDHG = [];
%%  PDHG.sigma  = '1';
%%  PDHG.tau    = '1.6e18';
%%  PDHG.theta  = '1';
%%  PDHG.lambda = '1e-3';
%%  PDHG.iter   = int2str(5);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  saveas(gcf, ['./figures/Example81_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);


%==============================
% S-PDHG
%==============================
%%  % S-PDHG **********************
%%  SPDHG = [];
%%  SPDHG.sigma  = '1';
%%  SPDHG.tau    = '2e17';
%%  SPDHG.theta  = '1';
%%  SPDHG.lambda = '1e-3';
%%  SPDHG.batch  = '100';
%%  SPDHG.epoch  = int2str(17);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection(pixelPressure, dx);
%%  a = axes;
%%  t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
%%  a.Visible = 'off'; 
%%  t.Visible = 'on'; 
%%  %saveas(gcf, ['./figures/Example81_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);


%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     DISTANCE ERROR                         ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
disp('******* PRIMAL DISTANCE ********');
disp('******* DUAL DISTANCE ********');
%======================================================================
% Gradient descent
%======================================================================
disp('GD');
GD.tau = {'5e17', '1e18', '2e18'};
GD.nIter = {30, 30, 30};
GD.lambda = '1e-4';
L = length(GD.tau);
clear GD_error_pd GD_error_dd;
for ii = 1:L
    disp(['STEP ', int2str(ii)])
    GD_error_pd{ii} = norm_distance(u0, 0*u0);
    GD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:GD.nIter{ii}
        disp(['iter', int2str(iter)])
        % Primal error
        ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pp)];
        % Dual error
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:GD.nIter{ii}, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 1 100]);
legend('GD 5e17', 'GD 1e18', 'GD 2e18');
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:GD.nIter{ii}, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 1e-19 1e-17])
legend('GD 5e17', 'GD 1e18', 'GD 2e18');
save ./results/error_vectors/GD_error_pd_lambda1em4 GD_error_pd GD;
save ./results/error_vectors/GD_error_dd_lambda1em4 GD_error_dd GD;
 
%======================================================================
% Stochastic Gradient descent
%======================================================================
disp('S-GD');
SGD.tau = {'2e18', '4e18', '8e18'};
SGD.nIter = {23, 30, 30};
SGD.lambda = '1e-4';
SGD.batch = '100';
L = length(SGD.tau);
clear SGD_error_pd SGD_error_dd;
for ii = 1:L
    disp(['STEP ', int2str(ii)])
    SGD_error_pd{ii} = norm_distance(u0, 0*u0);
    SGD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:SGD.nIter{ii}
        disp(['iter ', int2str(iter)])
        % Primal error
        ppmatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        SGD_error_pd{ii} = [SGD_error_pd{ii} norm_distance(u0, pp)];
        % Dual error
        tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SGD.nIter{ii}, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
axis([0 30 1 100]);
ax = gca;
ax.GridAlpha = 0.2;
legend('SGD 2e18', 'SGD 4e18', 'SGD 8e18');
% Plot dual
figure();
colors = winter(length(SGD.tau));
for ii = 1:length(SGD.tau)
    semilogy(0:SGD.nIter{ii}, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
axis([0 30 1e-19 1e-17]);
ax = gca;
ax.GridAlpha = 0.2;
legend('SGD 2e18', 'SGD 4e18', 'SGD 8e18');
save ./results/error_vectors/SGD_error_pd_lambda1em4 SGD_error_pd SGD;
save ./results/error_vectors/SGD_error_dd_lambda1em4 SGD_error_dd SGD;

%==============================
% FISTA
%==============================
disp('FISTA');
FISTA.tau = {'1e18'};
FISTA.nIter = {30};
FISTA.lambda = '1e-4';
L = length(FISTA.tau);
clear FISTA_error_pd FISTA_error_dd;
for ii = 1:L
    disp(ii)
    FISTA_error_pd{ii} = norm_distance(u0, 0*u0);
    FISTA_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:FISTA.nIter{ii}
        % Primal
        ppmatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        FISTA_error_pd{ii} = [FISTA_error_pd{ii} norm_distance(u0, pp)];
        % Dual
        tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        FISTA_error_dd{ii} = [FISTA_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot primal
figure();
colors = winter(length(FISTA.tau));
for ii = 1:length(FISTA.tau)
    semilogy(0:FISTA.nIter{ii}, FISTA_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
axis([0 30 1 100]);
ax = gca;
ax.GridAlpha = 0.2;
legend('FISTA 1e18');
% Plot dual
figure();
colors = winter(length(FISTA.tau));
for ii = 1:length(FISTA.tau)
    semilogy(0:FISTA.nIter{ii}, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
axis([0 30 1e-19 1e-17]);
ax = gca;
ax.GridAlpha = 0.2;
legend('FISTA 1e18');
save ./results/error_vectors/FISTA_error_pd_lambda1em4 FISTA_error_pd FISTA;
save ./results/error_vectors/FISTA_error_dd_lambda1em4 FISTA_error_dd FISTA;

%======================================================================
% PDHG
%======================================================================
disp('PDHG');
PDHG.tau = {'5e17', '1e18', '2e18'};
PDHG.nIter = {30, 30, 30};
PDHG.sigma = '1';
PDHG.theta = '1';
PDHG.lambda = '1e-4';
L = length(PDHG.tau);
clear PDHG_error_pd PDHG_error_dd; 
for ii = 1:L
    disp(['STEP ', int2str(ii)])
    PDHG_error_pd{ii} = norm_distance(u0, 0*u0);
    PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:PDHG.nIter{ii}
        disp(['iter ', int2str(iter)])
        % Primal
        ppmatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        PDHG_error_pd{ii} = [PDHG_error_pd{ii} norm_distance(u0, pp)];
        % Dual
        tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot Primal
figure();
colors = winter(length(PDHG.tau));
for ii = 1:length(PDHG.tau)
    semilogy(0:PDHG.nIter{ii}, PDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
axis([0 30 1 100]);
ax = gca;
ax.GridAlpha = 0.2;
legend('PDHG 5e17', 'PDHG 1e18', 'PDHG 2e18');
% Plot Dual
figure();
colors = winter(length(PDHG.tau));
for ii = 1:length(PDHG.tau)
    semilogy(0:PDHG.nIter{ii}, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
axis([0 30 1e-19 1e-17]);
ax = gca;
ax.GridAlpha = 0.2;
legend('PDHG 5e17', 'PDHG 1e18', 'PDHG 2e18');
save ./results/error_vectors/PDHG_error_pd_lambda1em4 PDHG_error_pd PDHG;
save ./results/error_vectors/PDHG_error_dd_lambda1em4 PDHG_error_dd PDHG;

%==============================
% SPDHG
%==============================
disp('S-PDHG');
SPDHG.tau = {'5e16', '1e17', '2e17'};
SPDHG.sigma = '1';
SPDHG.theta = '1';
SPDHG.batch = '100';
SPDHG.nIter = {30, 30, 30};
SPDHG.lambda = '1e-4';
L = length(SPDHG.tau);
clear SPDHG_error_pd SPDHG_error_dd;
for ii = 1:L
    disp(['STEP ', int2str(ii)])
    SPDHG_error_pd{ii} = norm_distance(u0, 0*u0);
    SPDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:SPDHG.nIter{ii}
        disp(['iter ', int2str(iter)])
        % Primal
        ppmatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        SPDHG_error_pd{ii} = [SPDHG_error_pd{ii} norm_distance(u0, pp)];
        % Dual
        tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SPDHG_error_dd{ii} = [SPDHG_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot
figure();
colors = winter(length(SPDHG.tau));
for ii = 1:length(SPDHG.tau)
    semilogy(0:SPDHG.nIter{ii}, SPDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 1 100]);
legend('SPDHG 5e16', 'SPDHG 1e16', 'SPDHG 2e17');
% Plot
figure();
colors = winter(length(SPDHG.tau));
for ii = 1:length(SPDHG.tau)
    semilogy(0:SPDHG.nIter{ii}, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 1e-19 1e-17]);
legend('SPDHG 5e16', 'SPDHG 1e16', 'SPDHG 2e17');
save ./results/error_vectors/SPDHG_error_pd_lambda1em4 SPDHG_error_pd SPDHG;
save ./results/error_vectors/SPDHG_error_dd_lambda1em4 SPDHG_error_dd SPDHG;

%======================================================================
% PLOT ALL PRIMAL
%======================================================================
load ./results/error_vectors/GD_error_pd_lambda1em4;
load ./results/error_vectors/SGD_error_pd_lambda1em4;
load ./results/error_vectors/FISTA_error_pd_lambda1em4;
load ./results/error_vectors/PDHG_error_pd_lambda1em4;
load ./results/error_vectors/SPDHG_error_pd_lambda1em4;
% Plot
figure();
semilogy(0:GD.nIter{2}, GD_error_pd{2}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{2}, SGD_error_pd{2}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{1}, FISTA_error_pd{1}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{2}, PDHG_error_pd{2}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{2}, SPDHG_error_pd{2}, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Primal Distance Error - homogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 1 100]);
saveas(gcf, ['./figures/Example81_pd_lambda1em4_error.fig']);

%======================================================================
% PLOT ALL DUAL
%======================================================================
load ./results/error_vectors/GD_error_dd_lambda1em4;
load ./results/error_vectors/SGD_error_dd_lambda1em4;
load ./results/error_vectors/FISTA_error_dd_lambda1em4;
load ./results/error_vectors/PDHG_error_dd_lambda1em4;
load ./results/error_vectors/SPDHG_error_dd_lambda1em4;

% Plot
figure();
semilogy(0:GD.nIter{2}, GD_error_dd{2}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{2}, SGD_error_dd{2}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{1}, FISTA_error_dd{1}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{2}, PDHG_error_dd{2}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{2}, SPDHG_error_dd{2}, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Dual Distance Error - homogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 1e-19 1e-17]);
saveas(gcf, ['./figures/Example81_dd_error.fig']);


%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                       LONG RECONSTRUCTION                  ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
disp('******* LONG RECONSTRUCTION ****');
GD.tau = {'4e17', '8e17'};
GD.lambda = '1e-3';
nIter = 100;
% PLOT ************************
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{2}, '_lambda', GD.lambda, '_iter', int2str(nIter), '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FB - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%******************************
L = length(GD.tau);
clear GD_error_dd_long;
for ii = 1:L
    disp(ii)
    GD_error_dd_long{ii} = norm_distance(y0, 0*y0);
    GD_error_pd_long{ii} = norm_distance(u0, 0*u0);
    for iter = 1:100
        disp(['iter ', iter])
        % Primal error
        ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        GD_error_pd_long{ii} = [GD_error_pd_long{ii} norm_distance(u0, pp)];
        % Dual error
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_dd_long{ii} = [GD_error_dd_long{ii} norm_distance(y0, yi)];
    end
end

load ./results/error_vectors/GD_error_long_reconstruction;
% Primal Error plot
figure();
colors = winter(length(GD.tau));
for ii = 1:length(GD.tau)
    semilogy(0:100, GD_error_pd_long{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
grid on;
ax = gca;
ax.GridAlpha = 0.2;
box on;
% Dual Error Plot
figure();
colors = winter(length(GD.tau));
for ii = 1:length(GD.tau)
    semilogy(0:100, GD_error_dd_long{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
grid on;
ax = gca;
ax.GridAlpha = 0.2;
box on;

% Save results
%save results/error_vectors/GD_error_long_reconstruction GD_error_dd_long GD_error_pd_long GD;
