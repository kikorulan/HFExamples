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

%============================================================================================================================================
% DATA
%============================================================================================================================================
%==================================================
% INITIAL PRESSURE
%==================================================
% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);

%==================================================
% TIME SIGNAL - FORWARD DATA
%==================================================
% Import data
filenameData = './input_data/forwardSignal_reference_14400sensors.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRD = timeSignal(1, :);
y0 = timeSignal(2:end, :);

%========================================================================================================================
% PRIMAL ERROR
%========================================================================================================================
nIter = 5;
%====================
% Gradient descent - lambda variation
%====================
disp('GD');
GD.tau    = {'1e18'};
GD.lambda = {'1e-2'};
clear GD_error_pd;
for ii = 1:length(GD.lambda)
    disp(ii)
    GD_error_pd{ii} = norm_distance(u0, 0*u0);
    for iter = 1:nIter
        pixelPressure = importdata(['./results/adjoint/FB/pixelPressure_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
        GD_error_pd{ii} = [GD_error_pd{ii} norm_distance(u0, pixelPressure)];
    end
end

%====================
% Stochastic Gradient descent - lambda variation
%====================
disp('S-GD');
%SGD.batch  = {'400', '800'};
SGD.batch  = {'90', '200', '1600', '3200', '6400'};
SGD.tau    = {'8e18'};
SGD.lambda = {'1e-4', '2e-4', '5e-4', '1e-3', '2e-3', '5e-3'};
clear SGD_error_pd;
for ibatch = 1:length(SGD.batch)
    for ii = 1:length(SGD.lambda)
        disp([ibatch, ii])
        SGD_error_pd{ibatch}{ii} = norm_distance(u0, 0*u0);
        for iter = 1:nIter
            pixelPressure = importdata(['./results/adjoint/S-FB/', SGD.batch{ibatch}, '/pixelPressure_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch{ibatch}, '_subepoch', int2str(iter), '.dat'], ' ', 0);
            pixelPressure = max(0, matrix2cube(pixelPressure, Nz));
            SGD_error_pd{ibatch}{ii} = [SGD_error_pd{ibatch}{ii} norm_distance(u0, pixelPressure)];
        end
    end
end
% Plot
x_axis = 0:nIter;
colors = winter(length(SGD.lambda));
for ibatch = 1:length(SGD.batch)
    figure();
    semilogy(x_axis, GD_error_pd{1}, 'Color', 'r', 'Linewidth', 1.5);
    hold on;
    for ii = 1:length(SGD.lambda)
        semilogy(x_axis, SGD_error_pd{ibatch}{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    end
    axis([0 nIter 1 200]);
    legend('GD', 'S-GD lambda = 1e-4', 'S-GD lambda = 2e-4', 'S-GD lambda = 5e-4', 'S-GD lambda = 1e-3', 'S-GD lambda = 2e-3', 'S-GD lambda = 5e-3');
    title(['Primal Distance Error , batch = ', SGD.batch{ibatch}]);
    grid on;
    box on;
    ax = gca;
    ax.GridAlpha = 0.2;
end


%========================================================================================================================
% DUAL DISTANCE
%========================================================================================================================
nIter = 5;
%====================
% Gradient descent - lambda variation
%====================
disp('GD');
GD.tau    = {'1e18'};
GD.lambda = {'1e-2'};
clear GD_error_dd;
for ii = 1:length(GD.lambda)
    disp(ii)
    GD_error_dd{ii} = norm_distance(y0, 0*y0);
    for iter = 1:nIter
        forwardSignal = importdata(['./results/forward/FB/forwardSignal_GD_tau', GD.tau{1}, '_lambda', GD.lambda{ii}, '_iter', int2str(iter), '.dat'], ' ', 0);
        forwardSignal = forwardSignal(2:end, :);
        GD_error_dd{ii} = [GD_error_dd{ii} norm_distance(y0, forwardSignal)];
    end
end

%====================
% Stochastic Gradient descent - lambda variation
%====================
disp('S-GD');
SGD.batch  = {'90', '200'};
SGD.tau    = {'8e18'};
SGD.lambda = {'1e-4', '2e-4', '5e-4', '1e-3', '2e-3', '5e-3'};
clear SGD_error_dd;
for ibatch = 1:length(SGD.batch)
    for ii = 1:length(SGD.lambda)
        disp([ibatch, ii])
        SGD_error_dd{ibatch}{ii} = norm_distance(y0, 0*y0);
        for iter = 1:nIter
            forwardSignal = importdata(['./results/forward/S-FB/', SGD.batch{ibatch}, '/forwardSignal_S-GD_tau', SGD.tau{1}, '_lambda', SGD.lambda{ii}, '_batch', SGD.batch{ibatch}, '_subepoch', int2str(iter), '.dat'], ' ', 0);
            forwardSignal = forwardSignal(2:end, :);
            SGD_error_dd{ibatch}{ii} = [SGD_error_dd{ibatch}{ii} norm_distance(y0, forwardSignal)];
        end
    end
end
% Plot
x_axis = 0:nIter;
colors = winter(length(SGD.lambda));
for ibatch = 1:length(SGD.batch)
    figure();
    semilogy(x_axis, GD_error_dd{1}, 'Color', 'r', 'Linewidth', 1.5);
    hold on;
    for ii = 1:length(SGD.lambda)
        semilogy(x_axis, SGD_error_dd{ibatch}{ii}, 'Color', colors(ii, :), 'Linewidth', 1.5);
    end
    axis([0 nIter 1e-18 2e-17]);
    legend('GD', 'S-GD lambda = 1e-4', 'S-GD lambda = 2e-4', 'S-GD lambda = 5e-4', 'S-GD lambda = 1e-3', 'S-GD lambda = 2e-3', 'S-GD lambda = 5e-3');
    title(['Dual Distance Error , batch = ', SGD.batch{ibatch}]);
    grid on;
    box on;
    ax = gca;
    ax.GridAlpha = 0.2;
end
