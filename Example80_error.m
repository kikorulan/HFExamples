% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex80_3D_veins_subsampled;

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
%saveas(gcf, './figures/Example80_initial_pressure.fig');

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
% GD **************************
GD = [];
GD.tau    = '4e17';
GD.lambda = '1e-2';
GD.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['GD - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.fig']);

%==============================
% Stochastic Gradient Descent
%==============================
% S-GD ************************
SGD = [];
SGD.tau    = '4e18';
SGD.lambda = '3e-4';
SGD.batch  = '90';
SGD.epoch  = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', batch = ', SGD.batch, ', epoch = ', SGD.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', SGD.epoch, '.fig']);

%==============================
% FISTA
%==============================
% FISTA ***********************
FISTA = [];
FISTA.tau    = '4e17';
FISTA.lambda = '1e-2';
FISTA.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['FISTA - t = ', FISTA.tau, ', l = ', FISTA.lambda, ', iter = ', FISTA.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.fig']);

%==============================
% PDHG
%==============================
% PDHG ************************
PDHG = [];
PDHG.sigma  = '1';
PDHG.tau    = '8e17';
PDHG.theta  = '1';
PDHG.lambda = '1e-2';
PDHG.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['PDHG - s = ', PDHG.sigma, ', t = ', PDHG.tau, ', l = ', PDHG.lambda, ', iter = ', PDHG.iter, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.fig']);

%==============================
% S-PDHG
%==============================
% S-PDHG **********************
SPDHG = [];
SPDHG.sigma  = '1';
SPDHG.tau    = '8e18';
SPDHG.theta  = '1';
SPDHG.lambda = '1e-3';
SPDHG.batch  = '90';
SPDHG.epoch  = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['S-PDHG - s = ', SPDHG.sigma, ', t = ', SPDHG.tau, ', l = ', SPDHG.lambda, ', batch = ', SPDHG.batch, ', epoch = ', SPDHG.epoch, ' - homogeneous SS']);
a.Visible = 'off'; 
t.Visible = 'on'; 
%saveas(gcf, ['./figures/Example80_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', SPDHG.epoch, '.fig']);


%========================================================================================================================
% PRIMAL DISTANCE ERROR
%========================================================================================================================
disp('******* PRIMAL DISTANCE ********');
nIter = 5;
x_axis = 0:nIter;
%==============================
% Gradient descent
%==============================
disp('GD');
GD.tau = {'2e17', '4e17', '8e17', '1.6e18', '3.2e18', '6.4e18'};
L = length(GD.tau);
clear GD_error_pd;
for ii = 1:L
    disp(ii)
    GD_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:nIter
        ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pp)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax.GridAlpha = 0.2;
legend('GD 2e17', 'GD 4e17', 'GD 8e17', 'GD 1.6e18', 'GD 3.2e18', 'GD 6.4e18');

%==============================
% Stochastic Gradient descent
%==============================
disp('S-GD');
SGD.tau = {'2e18', '4e18', '8e18'};
L = length(SGD.tau);
for ii = 1:L
    disp(ii)
    SGD_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:nIter
        ppmatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        SGD_error_pd{ii} = [SGD_error_pd{ii} norm_distance(u0, pp)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax.GridAlpha = 0.2;
legend('SGD 2e18', 'SGD 4e18', 'SGD 8e18');

%==============================
% FISTA
%==============================
disp('FISTA');
FISTA.tau = {'2e17', '4e17', '8e17', '1.6e18', '3.2e18', '6.4e18'};
L = length(FISTA.tau);
clear FISTA_error_pd;
for ii = 1:L
    disp(ii)
    FISTA_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:nIter
        ppmatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        FISTA_error_pd{ii} = [FISTA_error_pd{ii} norm_distance(u0, pp)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, FISTA_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax.GridAlpha = 0.2;
legend('FISTA 2e17', 'FISTA 4e17', 'FISTA 8e17', 'FISTA 1.6e18', 'FISTA 3.2e18', 'FISTA 6.4e18');

%==============================
% PDHG
%==============================
disp('PDHG');
PDHG.tau = {'4e17', '8e17', '1.6e18', '3.2e18', '6.4e18'};
L = length(PDHG.tau);
clear PDHG_error_pd;
for ii = 1:L
    disp(ii)
    PDHG_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:5
        ppmatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        PDHG_error_pd{ii} = [PDHG_error_pd{ii} norm_distance(u0, pp)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, PDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax.GridAlpha = 0.2;
legend('PDHG 4e17', 'PDHG 8e17', 'PDHG 1.6e18', 'PDHG 3.2e18', 'PDHG 6.4e18');

%==============================
% SPDHG
%==============================
disp('S-PDHG');
SPDHG.sigma = '2';
SPDHG.tau = {'4e18', '8e18', '1.6e19'};
L = length(SPDHG.tau);
clear SPDHG_error_pd;
for ii = 1:L
    disp(ii)
    SPDHG_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:nIter
        ppmatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        SPDHG_error_pd{ii} = [SPDHG_error_pd{ii} norm_distance(u0, pp)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, SPDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 10 100]);
legend('SPDHG 4e18', 'SPDHG 8e18', 'SPDHG 1.6e19');


%==============================
% PLOT ALL
%==============================
% Plot
figure();
%semilogy(x_axis, GD_error_pd{4}, 'Color', 'r', 'Linewidth', 1.5);
hold on;
%semilogy(x_axis, SGD_error_pd{2}, 'Color', 'g', 'Linewidth', 1.5);
%semilogy(x_axis, FISTA_error_pd{4}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(x_axis, PDHG_error_pd{4}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(x_axis, SPDHG_error_pd{2}, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Primal Distance Error - homogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 nIter 10 100]);
%saveas(gcf, ['./figures/Example80_pd_error.fig']);

%========================================================================================================================
% DUAL DISTANCE ERROR
%========================================================================================================================
disp('******* DUAL DISTANCE ********');
%==============================
% Gradient descent
%==============================
disp('GD');
GD.tau = {'2e17', '4e17', '8e17', '1.6e18', '3.2e18', '6.4e18'};
L = length(GD.tau);
clear GD_error_dd;
for ii = 1:L
    disp(ii)
    GD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
legend('GD 2e17', 'GD 4e17', 'GD 8e17', 'GD 1.6e18', 'GD 3.2e18', 'GD 6.4e18');

%==============================
% Stochastic Gradient descent
%==============================
disp('S-GD');
SGD.tau = {'2e18', '4e18', '8e18'};
L = length(SGD.tau);
clear SGD_error_dd;
for ii = 1:L
    disp(ii)
    SGD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
legend('SGD 2e18', 'SGD 4e18', 'SGD 8e18');


%==============================
% FISTA
%==============================
disp('FISTA');
FISTA.tau = {'2e17', '4e17', '8e17', '1.6e18', '3.2e18', '6.4e18'};
L = length(FISTA.tau);
clear FISTA_error_dd;
for ii = 1:L
    disp(ii)
    FISTA_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        FISTA_error_dd{ii} = [FISTA_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
legend('FISTA 2e17', 'FISTA 4e17', 'FISTA 8e17', 'FISTA 1.6e18', 'FISTA 3.2e18', 'FISTA 6.4e18');


%==============================
% PDHG
%==============================
disp('PDHG');
PDHG.tau = {'4e17', '8e17', '1.6e18', '3.2e18', '6.4e18'};
L = length(PDHG.tau);
clear PDHG_error_dd;
for ii = 1:L
    disp(ii)
    PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        %tSignal = importdata(['./output_data/a2.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
legend('PDHG 4e17', 'PDHG 8e17', 'PDHG 1.6e18', 'PDHG 3.2e18', 'PDHG 6.4e18');


%==============================
% SPDHG
%==============================
disp('S-PDHG');
SPDHG.tau = {'4e18', '8e18', '1.6e19'};
L = length(SPDHG.tau);
clear SPDHG_error_dd;
for ii = 1:L
    disp(ii)
    SPDHG_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:5
        tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0); 
        yi = tSignal(2:end, :);
        SPDHG_error_dd{ii} = [SPDHG_error_dd{ii} norm_distance(y0, yi)];
    end
end
% Plot
figure();
colors = winter(L);
for ii = 1:L
    semilogy(x_axis, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
legend('SPDHG 4e18', 'SPDHG 8e18', 'SPDHG 1.6e19');


%==============================
% PLOT ALL
%==============================
% Plot
figure();
semilogy(x_axis, GD_error_dd{3}, 'Color', 'r', 'Linewidth', 1.5);
hold on;
semilogy(x_axis, SGD_error_dd{2}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(x_axis, FISTA_error_dd{3}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(x_axis, PDHG_error_dd{3}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(x_axis, SPDHG_error_dd{1}, 'Color', 'c', 'Linewidth', 1.5);
legend('GD', 'S-GD', 'FISTA', 'PDHG', 'S-PDHG');
%legend('GD', 'FISTA', 'PDHG', 'S-PDHG');
title('Dual Distance Error - homogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%saveas(gcf, ['./figures/Example80_dd_error.fig']);



%========================================================================================================================
% MEMORY VERIFICATION
%========================================================================================================================
disp('******* PRIMAL DISTANCE ********');
nIter = 5;
x_axis = 0:nIter;
%==============================
% Gradient descent
%==============================
disp('GD');
GD.tau = '1.6e18';
clear GD_error_pd;
GD_error_pd_master = norm_distance(u0, 0*u0);
GD_error_pd_repair = norm_distance(u0, 0*u0);
for iter = 1:nIter
    ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    GD_error_pd_master = [GD_error_pd_master norm_distance(u0, pp)];
    ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '_mem.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    GD_error_pd_repair = [GD_error_pd_repair norm_distance(u0, pp)];
end
% Plot
figure();
semilogy(x_axis, GD_error_pd_master, 'Color', 'b', 'Linewidth', 3)
hold on;
semilogy(x_axis, GD_error_pd_repair, 'Color', 'r', 'Linewidth', 1.5)
ax.GridAlpha = 0.2;
legend('master', 'repair');
