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
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%========================================================================================================================
% PRIMAL AND DUAL DATA
%========================================================================================================================
% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);

% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_14400sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 5;
%============================================================
% PARAMETERS
%============================================================
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
PDHG.sigma  = '2e1';
PDHG.tau    = '1e18';
PDHG.theta  = '1';
PDHG.lambda = '1e-2';
PDHG.iter   = int2str(iter);
% S-PDHG *********************************
SPDHG = [];
SPDHG.sigma  = '1';
SPDHG.tau    = '1.5e19';
SPDHG.theta  = '1';
SPDHG.lambda = '3e-4';
SPDHG.batch  = '90';
SPDHG.epoch  = int2str(iter);


%========================================================================================================================
% PRIMAL DISTANCE ERROR
%========================================================================================================================
disp('******* PRIMAL DISTANCE ********');
% Gradient descent
disp('GD');
GD_error_pd = norm_distance(u0, 0*u0);
for iter = 1:5
    ppmatrix = importdata(['./results/adjoint/GD/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    GD_error_pd = [GD_error_pd norm_distance(u0, pp)];
end

% Stochastic Gradient descent
disp('S-GD');
SGD_error_pd = norm_distance(u0, 0*u0);
for iter = 1:5
    ppmatrix = importdata(['./results/adjoint/S-GD/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    SGD_error_pd = [SGD_error_pd norm_distance(u0, pp)];
end

% FISTA
disp('FISTA');
FISTA_error_pd = norm_distance(u0, 0*u0);
for iter = 1:5
    ppmatrix = importdata(['./results/adjoint/FISTA/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
    pp = max(0, matrix2cube(ppmatrix, Nz));
    FISTA_error_pd = [FISTA_error_pd norm_distance(u0, pp)];
end

% PDHG
PDHG.length = 1;
PDHG.tau    = '2e18';
PDHG.sigma  = {'1'}; % {'1e-1', '2e-1', '5e-1', '1', '2', '5', '1e1', '2e1', '5e1'} ;
PDHG_error_pd = [];
for sim = 1:PDHG.length
    disp(['PDHG', int2str(sim)]);
    PDHG_error_pd{sim} = norm_distance(u0, 0*u0);
    for iter = 1:5
        ppmatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma{sim}, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        PDHG_error_pd{sim} = [PDHG_error_pd{sim} norm_distance(u0, pp)];
    end
end

% S-PDHG
SPDHG.length = 6;
SPDHG.tau    = '1.5e19';
SPDHG.sigma  = {'1e-1', '2e-1', '5e-1', '1', '2', '5', '1e1', '2e1', '5e1'};
SPDHG_error_pd = [];
for sim = 1:SPDHG.length
    disp(['SPDHG', int2str(sim)]);
    SPDHG_error_pd{sim} = norm_distance(u0, 0*u0);
    for iter = 1:5
        ppmatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma{sim}, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_epoch', int2str(iter), '.dat'], ' ', 0);
        pp = max(0, matrix2cube(ppmatrix, Nz));
        SPDHG_error_pd{sim} = [SPDHG_error_pd{sim} norm_distance(u0, pp)];
    end
end


% Plot
colorsPDHG = winter(PDHG.length);
colorsSPDHG = summer(SPDHG.length);
x_axis = [0 1 2 3 4 5];
figure();
semilogy(x_axis, GD_error_pd, 'Color', 'r', 'Linewidth', 1.5);
hold on;
semilogy(x_axis, SGD_error_pd, 'Color', 'g', 'Linewidth', 1.5);
semilogy(x_axis, FISTA_error_pd, 'Color', 'm', 'Linewidth', 1.5);
for sim = 1:PDHG.length
    semilogy(x_axis, PDHG_error_pd{sim}, 'Color', colorsPDHG(sim, :), 'Linewidth', 1.5);
end
for sim = 1:4
    semilogy(x_axis, SPDHG_error_pd{sim}, 'Color', colorsSPDHG(sim, :), 'Linewidth', 1.5);
end
       % ['PDHG - sigma = ', PDHG.sigma{2}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{3}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{4}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{5}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{6}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{7}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{8}, ', tau = ', PDHG.tau], ...
       % ['PDHG - sigma = ', PDHG.sigma{9}, ', tau = ', PDHG.tau], ...
legend('GD', 'S-GD', 'FISTA', ...
       ['PDHG - sigma = ', PDHG.sigma{1}, ', tau = ', PDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{1}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{2}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{3}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{4}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{5}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{6}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{7}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{8}, ', tau = ', SPDHG.tau], ...
       ['SPDHG - sigma = ', SPDHG.sigma{9}, ', tau = ', SPDHG.tau])
title('Primal Distance Error - heterogeneous SS');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%saveas(gcf, ['./figures/Example74_pd_error.fig']);
