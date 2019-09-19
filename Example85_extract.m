% Read data from files
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex85_3D_veins_subsampled;
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex85_3D_veins_subsampled;

clear all;
close all;

addpath('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_2DRT/VariationalMethods')
addpath('/cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/HighFreq_2DRT/utils')

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
   
%==================================================
% PRIMAL AND DUAL DATA
%==================================================
% Load Initial Pressure
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
  
% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_noisy5_3600sensors.dat'], ' ', 0);
y0 = time_signal(2:end, :);
  
%==========================================================================================================================================================================
%===============================                                                            ===============================================================================
%===============================                     DISTANCE ERROR                         ===============================================================================
%===============================                                                            ===============================================================================
%==========================================================================================================================================================================
disp('******* PRIMAL DISTANCE ********');
disp('******* DUAL DISTANCE ********');

%======================================================================
% PARAMETRIZATION
%======================================================================
% Gradient Descent
GD.extract = 0;
GD.tau = {'4e1', '8e1'};
GD.nIter = {26, 200};
GD.lambda = '1e-4';

% Stochastic Gradient Descent
SGD.extract = 0;
SGD.tau = {'8e1', '1.6e2', '3.2e2'};
SGD.nIter = {12, 30, 30};
SGD.lambda = '1e-4';
SGD.batch = '1800';

% FISTA
FISTA.extract = 1;
FISTA.tau = {'2e1', '4e1', '8e1'};
FISTA.nIter = {13, 40, 12};
FISTA.lambda = '1e-4';

% PDHG
PDHG.extract = 0;
PDHG.tau = {'1e1', '2e1', '4e1'};
PDHG.sigma = '1';
PDHG.nIter = {30, 30, 200};
PDHG.lambda = '1e-4';

% SPDHG
SPDHG.extract = 0;
SPDHG.tau = {'1e1', '2e1', '4e1'};
SPDHG.sigma = '2e-1';
SPDHG.nIter = {30, 30, 30};
SPDHG.lambda = '1e-4';
SPDHG.batch = '100';

%======================================================================
% Gradient descent
%======================================================================
if (GD.extract == 1)
disp('GD');
clear GD_error_pd GD_error_dd GD_error_data GD_error_reg;
for ii = 1:length(GD.tau)
    disp(ii)
    GD_error_pd{ii}   = norm_distance(u0, 0*u0);
    GD_error_data{ii} = obj_data(y0, 0*y0);
    GD_error_reg{ii}  = obj_reg(0, 0*u0); 
    GD_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:GD.nIter{ii}-1
        disp(['    iter ', int2str(iter)]) 
        % Primal error
        ppmatrix = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{ii}, '_lambda', GD.lambda, '_iter', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        GD_error_data{ii} = [GD_error_data{ii} obj_data(y0, yi)];
        GD_error_reg{ii}  = [GD_error_reg{ii} obj_reg(str2double(GD.lambda), pp)];
        GD_error_dd{ii}   = [GD_error_dd{ii} obj_function(y0, yi, str2double(GD.lambda), pp)];
    end
end
save ./results/error_vectors/GD_error_lambda1em4 GD_error_pd GD_error_data GD_error_reg GD_error_dd GD;
end

%======================================================================
% Stochastic Gradient descent
%======================================================================
if (SGD.extract == 1)
disp('SGD');
clear SGD_error_pd SGD_error_dd SGD_error_data SGD_error_reg;
for ii = 1:length(SGD.tau)
    disp(ii)
    SGD_error_pd{ii}   = norm_distance(u0, 0*u0);
    SGD_error_data{ii} = obj_data(0, 0*u0);
    SGD_error_reg{ii}  = obj_reg(0, 0*u0);
    SGD_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:SGD.nIter{ii}-1
        disp(['    iter ', int2str(iter)]) 
        % Primal error
        ppmatrix = importdata(['./results/adjoint/SFB/pixelPressure_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        SGD_error_pd{ii} = [SGD_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results/forward/SFB/forwardSignal_S-GD_tau', SGD.tau{ii}, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SGD_error_data{ii} = [SGD_error_data{ii} obj_data(y0, yi)];
        SGD_error_reg{ii}  = [SGD_error_reg{ii} obj_reg(str2double(SGD.lambda), pp)];
        SGD_error_dd{ii}   = [SGD_error_dd{ii} obj_function(y0, yi, str2double(SGD.lambda), pp)];
    end
end
save ./results/error_vectors/SGD_error_lambda1em4_batch1800 SGD_error_pd SGD_error_data SGD_error_reg SGD_error_dd SGD;
end

%======================================================================
% FISTA
%======================================================================
if (FISTA.extract == 1)
disp('FISTA');
clear FISTA_error_pd FISTA_error_dd FISTA_error_data FISTA_error_reg;
for ii = 1:length(FISTA.tau)
    disp(ii)
    FISTA_error_pd{ii}   = norm_distance(u0, 0*u0);
    FISTA_error_data{ii} = obj_data(y0, 0*y0);
    FISTA_error_reg{ii}  = obj_reg(0, 0*u0);
    FISTA_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:FISTA.nIter{ii}-1
        disp(['    iter ', int2str(iter)]) 
        % Primal error
        ppmatrix = importdata(['./results/adjoint/AFB/pixelPressure_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        FISTA_error_pd{ii} = [FISTA_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results/forward/AFB/forwardSignal_FISTA_tau', FISTA.tau{ii}, '_lambda', FISTA.lambda, '_iter', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        FISTA_error_data{ii} = [FISTA_error_data{ii} obj_data(y0, yi)];
        FISTA_error_reg{ii}  = [FISTA_error_reg{ii} obj_reg(str2double(FISTA.lambda), pp)];
        FISTA_error_dd{ii}   = [FISTA_error_dd{ii} obj_function(y0, yi, str2double(FISTA.lambda), pp)];
    end
end
save ./results/error_vectors/FISTA_error_lambda1em4 FISTA_error_pd FISTA_error_data FISTA_error_reg FISTA_error_dd FISTA;
end 

%======================================================================
% PDHG
%======================================================================
if (PDHG.extract == 1)
disp('PDHG');
clear PDHG_error_pd PDHG_error_dd PDHG_error_data PDHG_error_reg;
for ii = 1:length(PDHG.tau)
    disp(ii)
    PDHG_error_pd{ii}   = norm_distance(u0, 0*u0);
    PDHG_error_data{ii} = obj_data(y0, 0*y0);
    PDHG_error_reg{ii}  = obj_reg(0, 0*u0);
    PDHG_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:PDHG.nIter{ii}-1
        disp(['    iter ', int2str(iter)]) 
        % Primal error
        ppmatrix = importdata(['./results/adjoint/PDHG/pixelPressure_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta1_lambda', PDHG.lambda, '_iter', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        PDHG_error_pd{ii} = [PDHG_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results/forward/PDHG/forwardSignal_PDHG_sigma', PDHG.sigma, '_tau', PDHG.tau{ii}, '_theta1_lambda', PDHG.lambda, '_iter', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        PDHG_error_data{ii} = [PDHG_error_data{ii} obj_data(y0, yi)];
        PDHG_error_reg{ii}  = [PDHG_error_reg{ii} obj_reg(str2double(PDHG.lambda), pp)];
        PDHG_error_dd{ii}   = [PDHG_error_dd{ii} obj_function(y0, yi, str2double(PDHG.lambda), pp)];
    end
end
save ./results/error_vectors/PDHG_error_lambda1em4_sigma1 PDHG_error_pd PDHG_error_data PDHG_error_reg PDHG_error_dd PDHG;
end

%==============================
% SPDHG
%==============================
if (SPDHG.extract == 1)
disp('SPDHG');
clear SPDHG_error_pd SPDHG_error_dd SPDHG_error_data SPDHG_error_reg;
for ii = 1:length(SPDHG.tau)
    disp(ii)
    SPDHG_error_pd{ii}   = norm_distance(u0, 0*u0);
    SPDHG_error_data{ii} = obj_data(y0, 0*y0);
    SPDHG_error_reg{ii}  = obj_reg(0, 0*u0);
    SPDHG_error_dd{ii}   = obj_function(y0, 0*y0, 0, 0*u0);
    for iter = 1:SPDHG.nIter{ii}-1
        % Primal error
        ppmatrix = importdata(['./results/adjoint/SPDHG/pixelPressure_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta1_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter), '.dat'], ' ', 0);
        pp = matrix2cube(ppmatrix, Nz);
        pp_pos = max(0, pp);
        SPDHG_error_pd{ii} = [SPDHG_error_pd{ii} norm_distance(u0, pp_pos)];
        % Dual error
        tSignal = importdata(['./results/forward/SPDHG/forwardSignal_S-PDHG_sigma', SPDHG.sigma, '_tau', SPDHG.tau{ii}, '_theta1_lambda', SPDHG.lambda, '_batch', SPDHG.batch, '_subepoch', int2str(iter+1), '.dat'], ' ', 0);
        yi = tSignal(2:end, :);
        SPDHG_error_data{ii} = [SPDHG_error_data{ii} obj_data(y0, yi)];
        SPDHG_error_reg{ii}  = [SPDHG_error_reg{ii} obj_reg(str2double(SPDHG.lambda), pp)];
        SPDHG_error_dd{ii}   = [SPDHG_error_dd{ii} obj_function(y0, yi, str2double(SPDHG.lambda), pp)];
    end
end
save ./results/error_vectors/SPDHG_error_lambda1em4_sigma5em2_batch100 SPDHG_error_pd SPDHG_error_data SPDHG_error_reg SPDHG_error_dd SPDHG;
end

