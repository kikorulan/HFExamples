% Read data from files
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex86_3D_veins_het;
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex86_3D_veins_het;

clear all;
close all;

% Functions
[TV, D, DTV] = TVOperators(3, 'none');
norm_distance = @(x, y) sum((x(:) - y(:)).*(x(:) - y(:)));
obj_data = @(y0, y) 0.5*norm_distance(y0, y);
obj_reg  = @(lambda, u0) lambda*TV(u0);
obj_function = @(y0, y, lambda, u0) obj_data(y0, y) + obj_reg(lambda, u0);

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
plot_projection_compact(u0, dx);
%saveas(gcf, './figures/Example86_initial_pressure', 'epsc');
  
% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_noisy5_3600sensors.dat'], ' ', 0);
%time_signal = importdata(['./input_data/forwardSignal_reference_noisy5_3600sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);
figure;
imagesc(y0(1:60, :));
colorbar()
%saveas(gcf, './figures/Example86_forwardSignal_noisy5', 'epsc');
  
% Load Adjoint Pressure
%%  pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
%%  pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
%%  plot_projection_compact(pressure_adjoint, dx);
%%  %saveas(gcf, './figures/Example85_pixelPressure_adjoint', 'epsc');

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================

%==============================
% Gradient Descent
%==============================
% GD **************************
GD = [];
GD.tau    = '4e1';
GD.lambda = '1e-4';
GD.iter   = int2str(13);
%******************************
%pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressureMatrix = importdata(['./results/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx);
%saveas(gcf, ['./figures_paper/Example85_GD'], 'epsc');

%==============================
% Stochastic Gradient Descent
%==============================
%%  % S-GD ************************
%%  SGD = [];
%%  SGD.tau    = '1.6e2';
%%  SGD.lambda = '1e-4';
%%  SGD.batch  = '1800';
%%  SGD.epoch  = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/SFB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  %saveas(gcf, ['./figures_paper/Example80_S-GD'], 'epsc');

%==============================
% FISTA
%==============================
%%  % FISTA ***********************
%%  FISTA = [];
%%  FISTA.tau    = '4e1';
%%  FISTA.lambda = '1e-4';
%%  FISTA.iter   = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  %saveas(gcf, ['./figures_paper/Example80_FISTA'], 'epsc');

%==============================
% PDHG
%==============================
%%  % PDHG ************************
%%  PDHG = [];
%%  PDHG.sigma  = '0.5';
%%  PDHG.tau    = '1.6e2';
%%  PDHG.theta  = '1';
%%  PDHG.lambda = '1e-4';
%%  PDHG.iter   = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  %saveas(gcf, ['./figures_paper/Example80_PDHG'], 'epsc');

%==============================
% S-PDHG
%==============================
%%  % S-PDHG **********************
%%  SPDHG = [];
%%  SPDHG.sigma  = '1';
%%  SPDHG.tau    = '5';
%%  SPDHG.theta  = '1';
%%  SPDHG.lambda = '1e-4';
%%  SPDHG.batch  = '100';
%%  SPDHG.epoch  = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/SPDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  %saveas(gcf, ['./figures_paper/Example80_S-PDHG'], 'epsc');


%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     DISTANCE ERROR                         ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
disp('******* PRIMAL DISTANCE ********');
disp('******* DUAL DISTANCE ********');

% LOAD DATA
load ./results/error_vectors/GD_error_lambda1em4;
load ./results/error_vectors/SGD_error_lambda1em4_batch1800;
load ./results/error_vectors/FISTA_error_lambda1em4;
load ./results/error_vectors/PDHG_error_lambda1em4_sigma1;
load ./results/error_vectors/SPDHG_error_lambda1em4_sigma5em2_batch100;

%======================================================================
% Gradient descent
%======================================================================
disp('GD');
L = length(GD.tau);
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:GD.nIter{ii}-1, GD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
title('GD primal')
% Plot dual
figure();
for ii = 1:L
    semilogy(0:GD.nIter{ii}-1, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
axis([0 200 0.7 20]);
title('GD dual')
 
%======================================================================
% Stochastic Gradient descent
%======================================================================
disp('SGD');
L = length(SGD.tau);
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SGD.nIter{ii}-1, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 1 100]);
title('SGD primal')

% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SGD.nIter{ii}-1, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('SGD dual')
 

%======================================================================
% FISTA
%======================================================================
disp('FISTA');
L = length(FISTA.tau);
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:FISTA.nIter{ii}-1, FISTA_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 1 100]);
title('FISTA primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:FISTA.nIter{ii}-1, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('FISTA dual')
 

%======================================================================
% PDHG
%======================================================================
disp('PDHG');
L = length(PDHG.tau);
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:PDHG.nIter{ii}-1, PDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 1 100]);
title('PDHG primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:PDHG.nIter{ii}-1, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('PDHG dual')

%==============================
% SPDHG
%==============================
disp('SPDHG');
L = length(SPDHG.tau);
% Plot primal
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SPDHG.nIter{ii}-1, SPDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 1 100]);
title('SPDHG primal')
% Plot dual
figure();
colors = winter(L);
for ii = 1:L
    semilogy(0:SPDHG.nIter{ii}-1, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.7 20])
title('SPDHG dual')


%========================================================================================================================
% PLOT ALL
%========================================================================================================================
% Choose index
GD.plotIndex = 2;
SGD.plotIndex = 1;
FISTA.plotIndex = 1;
PDHG.plotIndex = 1;
SPDHG.plotIndex = 2;

%======================================================================
% PLOT PRIMAL
%======================================================================
% Plot
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_pd{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_pd{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_pd{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_pd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_pd{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 40 1e4]);
xlabel('iter/epoch');
ylabel('primal error');
set(gca,'FontSize',15);
%saveas(gcf, './figures_paper/Example86_pdError', 'epsc');

%======================================================================
% PLOT DUAL AND RELATIVE
%======================================================================
% Plot
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_dd{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_dd{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_dd{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_dd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_dd{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('dual error');
%set(gca,'FontSize',15);
%saveas(gcf, './figures_paper/Example80_ddError', 'epsc');

% Relative error
figure();
semilogy(1:GD.nIter{GD.plotIndex}-1, (GD_error_dd{GD.plotIndex}(1:end-1)-GD_error_dd{GD.plotIndex}(2:end))./GD_error_dd{GD.plotIndex}(2:end), 'Color', 'r', 'Linewidth', 1.5)
hold on;       
semilogy(1:SGD.nIter{SGD.plotIndex}-1, (SGD_error_dd{SGD.plotIndex}(1:end-1)-SGD_error_dd{SGD.plotIndex}(2:end))./SGD_error_dd{SGD.plotIndex}(2:end), 'Color', 'g', 'Linewidth', 1.5)
semilogy(1:FISTA.nIter{FISTA.plotIndex}-1, (FISTA_error_dd{FISTA.plotIndex}(1:end-1)-FISTA_error_dd{FISTA.plotIndex}(2:end))./FISTA_error_dd{FISTA.plotIndex}(2:end), 'Color', 'b', 'Linewidth', 1.5)
semilogy(1:PDHG.nIter{PDHG.plotIndex}-1, (PDHG_error_dd{PDHG.plotIndex}(1:end-1)-PDHG_error_dd{PDHG.plotIndex}(2:end))./PDHG_error_dd{PDHG.plotIndex}(2:end), 'Color', 'm', 'Linewidth', 1.5)
semilogy(1:SPDHG.nIter{SPDHG.plotIndex}-1, (SPDHG_error_dd{SPDHG.plotIndex}(1:end-1)-SPDHG_error_dd{SPDHG.plotIndex}(2:end))./SPDHG_error_dd{SPDHG.plotIndex}(2:end), 'Color', 'c', 'Linewidth', 1.5)
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([1 30 5e-3 1]);
xlabel('iter/epoch');
ylabel('relative error');
%set(gca,'FontSize',15);
%saveas(gcf, './figures_paper/Example80_relativeDualError', 'epsc');


%======================================================================
% PLOT REGULARIZATION AND DATA FIT
%======================================================================

% Plot Data Term
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_data{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_data{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_data{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_data{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_data{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('data term');

% Plot Regularization Term
figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, GD_error_reg{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, SGD_error_reg{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, FISTA_error_reg{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_reg{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, SPDHG_error_reg{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('regularization');


%======================================================================
% RELATIVE DISTANCE TO OBJECTIVE FUNCTION
%======================================================================
OBJ_VAL = 0.94610;
%OBJ_VAL = 0.93448;
rel_distance = @(curve, iter) (curve(iter) - OBJ_VAL)/(curve(1) - OBJ_VAL);

figure();
semilogy(0:GD.nIter{GD.plotIndex}-1, rel_distance(GD_error_dd{GD.plotIndex}, 1:GD.nIter{GD.plotIndex}), 'Color', 'r', 'Linewidth', 1.5);
hold on;       
semilogy(0:SGD.nIter{SGD.plotIndex}-1, rel_distance(SGD_error_dd{SGD.plotIndex}, 1:SGD.nIter{SGD.plotIndex}), 'Color', 'g', 'Linewidth', 1.5);
semilogy(0:FISTA.nIter{FISTA.plotIndex}-1, rel_distance(FISTA_error_dd{SGD.plotIndex}, 1:FISTA.nIter{FISTA.plotIndex}), 'Color', 'b', 'Linewidth', 1.5);
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, rel_distance(PDHG_error_dd{PDHG.plotIndex}, 1:PDHG.nIter{PDHG.plotIndex}), 'Color', 'm', 'Linewidth', 1.5);
semilogy(0:SPDHG.nIter{SPDHG.plotIndex}-1, rel_distance(SPDHG_error_dd{SPDHG.plotIndex}, 1:SPDHG.nIter{SPDHG.plotIndex}), 'Color', 'c', 'Linewidth', 1.5);
legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('data term');

title('Example 86 - Objective function')
axis([0 100 1e-6 1 ])
saveas(gcf, 'figures/Example86_objectiveFunction', 'png');



%%  % Plot as function of factor
%%  S_PDHG_func = @(factor) SPDHG_error_data{SPDHG.plotIndex}(31) + factor*SPDHG_error_reg{SPDHG.plotIndex}(31);
%%  GD_func = @(factor) GD_error_data{GD.plotIndex}(31) + factor*GD_error_reg{GD.plotIndex}(31);
%%  
%%  factor = [0 0.001 0.01 0.1 1 10 1000];
%%  figure;
%%  semilogy(S_PDHG_func(factor))
%%  hold on;
%%  semilogy(GD_func(factor))

