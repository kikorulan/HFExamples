% Read data from files
%cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled;
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex85_3D_veins_subsampled;

clear all;
close all;

% Functions
lambda = 1e-20;
[TV, D, DTV] = TVOperators(3, 'none');
norm_distance = @(x, y) sum((x(:) - y(:)).*(x(:) - y(:)));
obj_function = @(y0, y, u0) 0.5*norm_distance(y0, y) + lambda*TV(u0);

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

%========================================================================================================================
% ITERATIVE RECONSTRUCTION
%========================================================================================================================
iter = 30;

%==============================
% Gradient Descent
%==============================
% GD **************************
GD = [];
GD.tau    = '1.6e18';
GD.lambda = '1e-20';
GD.iter   = int2str(iter);
%******************************
pixelPressureMatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection_compact(pixelPressure, dx);
%saveas(gcf, ['./figures_paper/Example85_GD'], 'epsc');

%==============================
% Stochastic Gradient Descent
%==============================
%%  % S-GD ************************
%%  SGD = [];
%%  SGD.tau    = '3.2e18';
%%  SGD.lambda = '1e-4';
%%  SGD.batch  = '100';
%%  SGD.epoch  = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  saveas(gcf, ['./figures_paper/Example80_S-GD'], 'epsc');

%==============================
% FISTA
%==============================
%%  % FISTA ***********************
%%  FISTA = [];
%%  FISTA.tau    = '4e17';
%%  FISTA.lambda = '1e-4';
%%  FISTA.iter   = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau, '_lambda', FISTA.lambda, '_iter', FISTA.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  saveas(gcf, ['./figures_paper/Example80_FISTA'], 'epsc');

%==============================
% PDHG
%==============================
%%  % PDHG ************************
%%  PDHG = [];
%%  PDHG.sigma  = '1';
%%  PDHG.tau    = '8e17';
%%  PDHG.theta  = '1';
%%  PDHG.lambda = '1e-4';
%%  PDHG.iter   = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', PDHG.iter, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  saveas(gcf, ['./figures_paper/Example80_PDHG'], 'epsc');

%==============================
% S-PDHG
%==============================
%%  % S-PDHG **********************
%%  SPDHG = [];
%%  SPDHG.sigma  = '1';
%%  SPDHG.tau    = '1e17';
%%  SPDHG.theta  = '1';
%%  SPDHG.lambda = '1e-4';
%%  SPDHG.batch  = '100';
%%  SPDHG.epoch  = int2str(30);
%%  %******************************
%%  pixelPressureMatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', SPDHG.epoch, '.dat'], ' ', 0);
%%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
%%  plot_projection_compact(pixelPressure, dx);
%%  saveas(gcf, ['./figures_paper/Example80_S-PDHG'], 'epsc');


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
GD.tau = {'1e17', '2e17', '4e17', '8e17', '1.6e18'};
GD.nIter = {30, 30, 30, 30, 30};
GD.lambda = '1e-20';
L = length(GD.tau);
clear GD_error_pd GD_error_dd;
for ii = 1:L
    disp(ii)
    GD_error_pd{ii} = norm_distance(u0, 0*u0);
    GD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:GD.nIter{ii}
        % Primal error
        ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_dd{ii} = [GD_error_dd{ii} obj_function(y0, yi, pp)];
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
%axis([0 30 1 100]);
legend('GD 8e17');
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
%axis([0 30 1e-19 1e-17])
legend('GD 8e17');
save ./results/error_vectors/GD_error_pd_lambda1em20 GD_error_pd GD;
save ./results/error_vectors/GD_error_dd_lambda1em20 GD_error_dd GD;
 
%======================================================================
% Stochastic Gradient descent
%======================================================================
%%  disp('S-GD');
%%  SGD.tau = {'3.2e18'};
%%  SGD.nIter = {30};
%%  SGD.lambda = '1e-4';
%%  SGD.batch = '100';
%%  L = length(SGD.tau);
%%  clear SGD_error_pd SGD_error_dd;
%%  for ii = 1:L
%%      disp(ii)
%%      SGD_error_pd{ii} = norm_distance(u0, 0*u0);
%%      SGD_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:SGD.nIter{ii}
%%          % Primal error
%%          ppmatrix = importdata(['./results/adjoint/S-FB/pixelPressure_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          pp = max(0, matrix2cube(ppmatrix, Nz));
%%          SGD_error_pd{ii} = [SGD_error_pd{ii} norm_distance(u0, pp)];
%%          % Dual error
%%          tSignal = importdata(['./results/forward/S-FB/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          SGD_error_dd{ii} = [SGD_error_dd{ii} norm_distance(y0, yi)];
%%      end
%%  end
%%  % Plot primal
%%  figure();
%%  colors = winter(L);
%%  for ii = 1:L
%%      semilogy(0:SGD.nIter{ii}, SGD_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1 100]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('SGD 3.2e18');
%%  % Plot dual
%%  figure();
%%  colors = winter(length(SGD.tau));
%%  for ii = 1:length(SGD.tau)
%%      semilogy(0:SGD.nIter{ii}, SGD_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1e-19 1e-17]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('SGD 3.2e18');
%%  save ./results/error_vectors/SGD_error_pd_lambda1em4 SGD_error_pd SGD;
%%  save ./results/error_vectors/SGD_error_dd_lambda1em4 SGD_error_dd SGD;

%==============================
% FISTA
%==============================
%%  disp('FISTA');
%%  FISTA.tau = {'4e17'};
%%  FISTA.nIter = {30};
%%  FISTA.lambda = '1e-4';
%%  L = length(FISTA.tau);
%%  clear FISTA_error_pd FISTA_error_dd;
%%  for ii = 1:L
%%      disp(ii)
%%      FISTA_error_pd{ii} = norm_distance(u0, 0*u0);
%%      FISTA_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:FISTA.nIter{ii}
%%          % Primal
%%          ppmatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          pp = max(0, matrix2cube(ppmatrix, Nz));
%%          FISTA_error_pd{ii} = [FISTA_error_pd{ii} norm_distance(u0, pp)];
%%          % Dual
%%          tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          FISTA_error_dd{ii} = [FISTA_error_dd{ii} norm_distance(y0, yi)];
%%      end
%%  end
%%  % Plot primal
%%  figure();
%%  colors = winter(length(FISTA.tau));
%%  for ii = 1:length(FISTA.tau)
%%      semilogy(0:FISTA.nIter{ii}, FISTA_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1 100]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('FISTA 4e17');
%%  % Plot dual
%%  figure();
%%  colors = winter(length(FISTA.tau));
%%  for ii = 1:length(FISTA.tau)
%%      semilogy(0:FISTA.nIter{ii}, FISTA_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1e-19 1e-17]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('FISTA 4e17');
%%  save ./results/error_vectors/FISTA_error_pd_lambda1em4 FISTA_error_pd FISTA;
%%  save ./results/error_vectors/FISTA_error_dd_lambda1em4 FISTA_error_dd FISTA;

%======================================================================
% PDHG
%======================================================================
%%  disp('PDHG');
%%  PDHG.tau = {'8e17'};
%%  PDHG.nIter = {30};
%%  PDHG.sigma = '1';
%%  PDHG.theta = '1';
%%  PDHG.lambda = '1e-4';
%%  L = length(PDHG.tau);
%%  clear PDHG_error_pd PDHG_error_dd; 
%%  for ii = 1:L
%%      disp(ii)
%%      PDHG_error_pd{ii} = norm_distance(u0, 0*u0);
%%      PDHG_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:PDHG.nIter{ii}
%%          % Primal
%%          ppmatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          pp = max(0, matrix2cube(ppmatrix, Nz));
%%          PDHG_error_pd{ii} = [PDHG_error_pd{ii} norm_distance(u0, pp)];
%%          % Dual
%%          tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta', PDHG.theta, '_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          PDHG_error_dd{ii} = [PDHG_error_dd{ii} norm_distance(y0, yi)];
%%      end
%%  end
%%  % Plot Primal
%%  figure();
%%  colors = winter(length(PDHG.tau));
%%  for ii = 1:length(PDHG.tau)
%%      semilogy(0:PDHG.nIter{ii}, PDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 1 100]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('PDHG 8e17');
%%  % Plot Dual
%%  figure();
%%  colors = winter(length(PDHG.tau));
%%  for ii = 1:length(PDHG.tau)
%%      semilogy(0:PDHG.nIter{ii}, PDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  axis([0 30 5e-19 1e-17]);
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  legend('PDHG 8e17');
%%  save ./results/error_vectors/PDHG_error_pd_lambda1em4 PDHG_error_pd PDHG;
%%  save ./results/error_vectors/PDHG_error_dd_lambda1em4 PDHG_error_dd PDHG;

%==============================
% SPDHG
%==============================
%%  disp('S-PDHG');
%%  SPDHG.tau = {'5e16', '1e17', '2e17'};
%%  SPDHG.sigma = '1';
%%  SPDHG.theta = '1';
%%  SPDHG.batch = '100';
%%  SPDHG.nIter = {30, 30, 30};
%%  SPDHG.lambda = '1e-4';
%%  L = length(SPDHG.tau);
%%  clear SPDHG_error_pd SPDHG_error_dd;
%%  for ii = 1:L
%%      disp(ii)
%%      SPDHG_error_pd{ii} = norm_distance(u0, 0*u0);
%%      SPDHG_error_dd{ii} = norm_distance(y0, 0*y0);
%%      for iter = 1:SPDHG.nIter{ii}
%%          % Primal
%%          ppmatrix = importdata(['./results/adjoint/S-PDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          pp = max(0, matrix2cube(ppmatrix, Nz));
%%          SPDHG_error_pd{ii} = [SPDHG_error_pd{ii} norm_distance(u0, pp)];
%%          % Dual
%%          tSignal = importdata(['./results/forward/S-PDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta', SPDHG.theta, '_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
%%          yi = tSignal(2:end, :);
%%          SPDHG_error_dd{ii} = [SPDHG_error_dd{ii} norm_distance(y0, yi)];
%%      end
%%  end
%%  % Plot
%%  figure();
%%  colors = winter(length(SPDHG.tau));
%%  for ii = 1:length(SPDHG.tau)
%%      semilogy(0:SPDHG.nIter{ii}, SPDHG_error_pd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  axis([0 30 1 100]);
%%  legend('SPDHG 5e16', 'SPDHG 1e17', 'SPDHG 2e17');
%%  % Plot
%%  figure();
%%  colors = winter(length(SPDHG.tau));
%%  for ii = 1:length(SPDHG.tau)
%%      semilogy(0:SPDHG.nIter{ii}, SPDHG_error_dd{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5)
%%      hold on;
%%  end
%%  box on;
%%  grid on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  axis([0 30 1e-19 1e-17]);
%%  legend('SPDHG 5e16', 'SPDHG 1e17', 'SPDHG 2e17');
%%  save ./results/error_vectors/SPDHG_error_pd_lambda1em4 SPDHG_error_pd SPDHG;
%%  save ./results/error_vectors/SPDHG_error_dd_lambda1em4 SPDHG_error_dd SPDHG;

%======================================================================
% PLOT ALL PRIMAL
%======================================================================
%%  load ./results/error_vectors/GD_error_pd_lambda1em4;
%%  load ./results/error_vectors/SGD_error_pd_lambda1em4;
%%  load ./results/error_vectors/FISTA_error_pd_lambda1em4;
%%  load ./results/error_vectors/PDHG_error_pd_lambda1em4;
%%  load ./results/error_vectors/SPDHG_error_pd_lambda1em4;
%%  % Plot
%%  figure();
%%  semilogy(0:GD.nIter{1}, GD_error_pd{1}, 'Color', 'r', 'Linewidth', 1.5);
%%  hold on;       
%%  semilogy(0:SGD.nIter{1}, SGD_error_pd{1}, 'Color', 'g', 'Linewidth', 1.5);
%%  semilogy(0:FISTA.nIter{1}, FISTA_error_pd{1}, 'Color', 'b', 'Linewidth', 1.5);
%%  semilogy(0:PDHG.nIter{1}, PDHG_error_pd{1}, 'Color', 'm', 'Linewidth', 1.5);
%%  semilogy(0:SPDHG.nIter{2}, SPDHG_error_pd{2}, 'Color', 'c', 'Linewidth', 1.5);
%%  legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  axis([0 30 5 100]);
%%  xlabel('iter/epoch');
%%  ylabel('primal error');
%%  set(gca,'FontSize',15);
%%  saveas(gcf, './figures_paper/Example80_pdError', 'epsc');

%======================================================================
% PLOT ALL DUAL
%======================================================================
%%  load ./results/error_vectors/GD_error_dd_lambda1em4;
%%  load ./results/error_vectors/SGD_error_dd_lambda1em4;
%%  load ./results/error_vectors/FISTA_error_dd_lambda1em4;
%%  load ./results/error_vectors/PDHG_error_dd_lambda1em4;
%%  load ./results/error_vectors/SPDHG_error_dd_lambda1em4;


%%  % Choose index
%%  GD.plotIndex = 1;
%%  SGD.plotIndex = 1;
%%  FISTA.plotIndex = 1;
%%  PDHG.plotIndex = 1;
%%  SPDHG.plotIndex = 2;
%%  
%%  % Relative error
%%  figure();
%%  semilogy(1:GD.nIter{GD.plotIndex}, (GD_error_dd{GD.plotIndex}(1:end-1)-GD_error_dd{GD.plotIndex}(2:end))./GD_error_dd{GD.plotIndex}(2:end), 'Color', 'r', 'Linewidth', 1.5)
%%  hold on;       
%%  semilogy(1:SGD.nIter{SGD.plotIndex}, (SGD_error_dd{SGD.plotIndex}(1:end-1)-SGD_error_dd{SGD.plotIndex}(2:end))./SGD_error_dd{SGD.plotIndex}(2:end), 'Color', 'g', 'Linewidth', 1.5)
%%  semilogy(1:FISTA.nIter{FISTA.plotIndex}, (FISTA_error_dd{FISTA.plotIndex}(1:end-1)-FISTA_error_dd{FISTA.plotIndex}(2:end))./FISTA_error_dd{FISTA.plotIndex}(2:end), 'Color', 'b', 'Linewidth', 1.5)
%%  semilogy(1:PDHG.nIter{PDHG.plotIndex}, (PDHG_error_dd{PDHG.plotIndex}(1:end-1)-PDHG_error_dd{PDHG.plotIndex}(2:end))./PDHG_error_dd{PDHG.plotIndex}(2:end), 'Color', 'm', 'Linewidth', 1.5)
%%  semilogy(1:SPDHG.nIter{SPDHG.plotIndex}, (SPDHG_error_dd{SPDHG.plotIndex}(1:end-1)-SPDHG_error_dd{SPDHG.plotIndex}(2:end))./SPDHG_error_dd{SPDHG.plotIndex}(2:end), 'Color', 'c', 'Linewidth', 1.5)
%%  legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  axis([1 30 5e-3 1]);
%%  xlabel('iter/epoch');
%%  ylabel('relative error');
%%  set(gca,'FontSize',15);
%%  saveas(gcf, './figures_paper/Example80_relativeDualError', 'epsc');
%%  
%%  
%%  % Plot
%%  figure();
%%  semilogy(0:GD.nIter{GD.plotIndex}, GD_error_dd{GD.plotIndex}, 'Color', 'r', 'Linewidth', 1.5);
%%  hold on;       
%%  semilogy(0:SGD.nIter{SGD.plotIndex}, SGD_error_dd{SGD.plotIndex}, 'Color', 'g', 'Linewidth', 1.5);
%%  semilogy(0:FISTA.nIter{FISTA.plotIndex}, FISTA_error_dd{FISTA.plotIndex}, 'Color', 'b', 'Linewidth', 1.5);
%%  semilogy(0:PDHG.nIter{PDHG.plotIndex}, PDHG_error_dd{PDHG.plotIndex}, 'Color', 'm', 'Linewidth', 1.5);
%%  semilogy(0:SPDHG.nIter{SPDHG.plotIndex}, SPDHG_error_dd{SPDHG.plotIndex}, 'Color', 'c', 'Linewidth', 1.5);
%%  legend('FB', 'S-FB', 'AFB', 'PDHG', 'S-PDHG');
%%  grid on;
%%  box on;
%%  ax = gca;
%%  ax.GridAlpha = 0.2;
%%  axis([0 30 1e-19 1e-17]);
%%  xlabel('iter/epoch');
%%  ylabel('dual error');
%%  set(gca,'FontSize',15);
%%  saveas(gcf, './figures_paper/Example80_ddError', 'epsc');

