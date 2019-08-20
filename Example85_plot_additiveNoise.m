% Read data from files
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled;
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex85_3D_veins_subsampled;

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
%saveas(gcf, './figures/Example85_initial_pressure', 'epsc');
  
% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_noisy5_3600sensors.dat'], ' ', 0);
%time_signal = importdata(['./input_data/forwardSignal_reference_noisy5_3600sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);
figure;
imagesc(y0(1:60, :));
colorbar()
%saveas(gcf, './figures/Example85_forwardSignal_noisy5', 'epsc');
  
% Load Adjoint Pressure
pressure_adjoint = importdata('./input_data/pixelPressure_adjoint_3600sensors.dat', ' ', 0);
pressure_adjoint = matrix2cube(pressure_adjoint, Nz);
plot_projection_compact(pressure_adjoint, dx);
%saveas(gcf, './figures/Example85_pixelPressure_adjoint', 'epsc');

%%  % kWave
%%  p0_recon_PML = h5read('./output_data/Example85_adjoint_output_3600sensors.h5', '/p_final');
%%  PML_size = 10;
%%  pixelPressure_KW = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
%%  plot_projection_compact(pixelPressure_KW, 1);

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 30;

%==============================
% Gradient Descent
%==============================
%%  % GD **************************
%%  GD = [];
%%  GD.tau    = '8e1';
%%  GD.lambda = '1e-4';
%%  GD.iter   = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  %saveas(gcf, ['./figures_paper/Example85_GD'], 'epsc');

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
%%  %%  saveas(gcf, ['./figures_paper/Example80_S-GD'], 'epsc');

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

%======================================================================
% Gradient descent
%======================================================================
disp('GD');
GD.tau = {'4e1', '8e1', '1.6e2'};
GD.nIter = {30, 30, 30};
GD.lambda = '1e-4';
L = length(GD.tau);
clear GD_error_pd GD_error_dd GD_error_data GD_error_reg;
for ii = 1:L
    disp(ii)
    GD_error_pd{ii}   = norm_distance(u0, 0*u0);
    GD_error_data{ii} = obj_data(y0, 0*y0);
    GD_error_reg{ii}  = obj_reg(0, 0*u0); 
    GD_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:GD.nIter{ii}-1
        % Primal error
        ppmatrix = importdata(['./results_add/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results_add/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_data{ii} = [GD_error_data{ii} obj_data(y0, yi)];
        GD_error_reg{ii}  = [GD_error_reg{ii} obj_reg(str2double(GD.lambda), pp)];
        GD_error_dd{ii}   = [GD_error_dd{ii} obj_function(y0, yi, str2double(GD.lambda), pp)];
    end
end
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
%axis([0 30 5 20]);
legend('GD 8e17');
% Plot dual
figure();
for ii = 1:L
    semilogy(0:GD.nIter{ii}-1, GD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
    hold on;
end
box on;
grid on;
ax = gca;
axis([0 50 0.7 20]);
legend('GD 4e1', 'GD 8e1', 'GD 1.6e2');
save ./results_add/error_vectors/GD_error_lambda1em4 GD_error_pd GD_error_data GD_error_reg GD_error_dd GD;
 

%======================================================================
% PDHG
%======================================================================
disp('PDHG');
PDHG.tau = {'4e1', '8e1', '1.6e2'};
PDHG.sigma = '1';
PDHG.nIter = {30, 30, 30};
PDHG.lambda = '1e-4';
L = length(PDHG.tau);
clear PDHG_error_pd PDHG_error_dd PDHG_error_data PDHG_error_reg;
for ii = 1:L
    disp(ii)
    PDHG_error_pd{ii}   = norm_distance(u0, 0*u0);
    PDHG_error_data{ii} = obj_data(y0, 0*y0);
    PDHG_error_reg{ii}  = obj_reg(0, 0*u0);
    PDHG_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:PDHG.nIter{ii}-1
        % Primal error
        ppmatrix = importdata(['./results_add/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta1_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        PDHG_error_pd{ii} = [PDHG_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results_add/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta1_lambda', PDHG.lambda, '_iter', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        PDHG_error_data{ii} = [PDHG_error_data{ii} obj_data(y0, yi)];
        PDHG_error_reg{ii}  = [PDHG_error_reg{ii} obj_reg(str2double(PDHG.lambda), pp)];
        PDHG_error_dd{ii}   = [PDHG_error_dd{ii} obj_function(y0, yi, str2double(PDHG.lambda), pp)];
    end
end
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
%legend('PDHG 8e17');
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
save ./results_add/error_vectors/PDHG_error_lambda1em4_sigma1 PDHG_error_pd PDHG_error_data PDHG_error_reg PDHG_error_dd PDHG;

%========================================================================================================================
% PLOT ALL
%========================================================================================================================
load ./results_add/error_vectors/GD_error_lambda1em4;
load ./results_add/error_vectors/PDHG_error_lambda1em4_sigma1;

% Choose index
GD.plotIndex = 2;
PDHG.plotIndex = 2;

%======================================================================
% PLOT PRIMAL
%======================================================================
% Plot
figure();
loglog(0:GD.nIter{GD.plotIndex}-1, GD_error_pd{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
hold on;       
loglog(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_pd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
legend('FB', 'PDHG');
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
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_dd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
legend('FB', 'PDHG');
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
semilogy(1:PDHG.nIter{PDHG.plotIndex}-1, (PDHG_error_dd{PDHG.plotIndex}(1:end-1)-PDHG_error_dd{PDHG.plotIndex}(2:end))./PDHG_error_dd{PDHG.plotIndex}(2:end), 'Color', 'm', 'Linewidth', 1.5)
legend('FB', 'PDHG');
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
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_data{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
legend('FB', 'PDHG');
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
semilogy(0:PDHG.nIter{PDHG.plotIndex}-1, PDHG_error_reg{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
legend('FB', 'PDHG');
grid on;
box on;
ax = gca;
ax.GridAlpha = 0.2;
%axis([0 30 0.9 20]);
xlabel('iter/epoch');
ylabel('regularization');



%%  % Plot as function of factor
%%  S_PDHG_func = @(factor) SPDHG_error_data{SPDHG.plotIndex}(31) + factor*SPDHG_error_reg{SPDHG.plotIndex}(31);
%%  GD_func = @(factor) GD_error_data{GD.plotIndex}(31) + factor*GD_error_reg{GD.plotIndex}(31);
%%  
%%  factor = [0 0.001 0.01 0.1 1 10 1000];
%%  figure;
%%  semilogy(S_PDHG_func(factor))
%%  hold on;
%%  semilogy(GD_func(factor))

