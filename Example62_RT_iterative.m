%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex62_recon2D_subsample;


close all;
% Measure computational time
tic;
start_time = clock;

%load sensor_data.mat;
%load recon_data_sensors.mat;
load recon_data_adjoint.mat;
%load recon_data_RT.mat;

%==================================================================================
% SIMULATION
%==================================================================================
% Initialise forward and inverse data
if (exist('y0')); clear y0; end;
for n = 1:length(source)
    if(isempty(source(n).aForward_initial))
        source(n).setForwardSignal_initial(source(n).aForward);
    end
    y0(n, :) = source(n).aForward_initial;
end

%==============================
% FISTA
%==============================
x0 = Rgrid.pixelAReverse/max(Rgrid.pixelAReverse(:));
% Define parameters - FISTA
lambda_fista = 1e-15;
lipschitz = 2e-12;
nIter = 50;

% Full operator
forward_operator_full = @(initial_pressure) Rgrid.operator_forward(source, 'all', initial_pressure, norm_factor);
inverse_operator_full = @(forward_data) Rgrid.operator_inverse(source, 'all', forward_data);

% Run
pixelRT = fista(y0, x0, lambda_fista, lipschitz, nIter, forward_operator_full, inverse_operator_full);
%pixelRT = gradient_descent(b, norm_factor*x0, lambda, nIter, forward_operator, inverse_operator);

%%  %==============================
%%  % SPDHG
%%  %==============================
%%  x0 = Rgrid.pixelAReverse;
%%  % Define parameters - SPDHG
%%  lambda = 1;
%%  sigma = 1e5; %1e5
%%  tau = 1e11;
%%  theta = 1;
%%  nEpochs = 10;
%%  nSensors = 46;
%%  
%%  % Define operators
%%  %norm_factor = max(Rgrid.u(:))/max(Rgrid.pixelAReverse(:));
%%  norm_factor = 1;
%%  forward_operator = @(index, initial_pressure) Rgrid.operator_forward(source, index, initial_pressure, norm_factor);
%%  inverse_operator = @(index, forward_data) Rgrid.operator_inverse(source, index, forward_data);
%%  
%%  % Parameters
%%  para.maxIter = 20;
%%  para.constraint = 'positivity';
%%  
%%  % Proximal operators
%%  proxG = @(x, a) conTVdenoising(x, 5e-3, para);
%%  proxF = @(y, a, b) (y-a*y0(b, :))/(1+a); 
%%  
%%  % Run 
%%  pixelRT_SPDHG = StochasticPDHG(x0, y0, proxF, proxG, forward_operator, inverse_operator, sigma, tau, theta, nEpochs, nSensors);

%%  %==============================
%%  % k-Wave
%%  %==============================
%%  % Define parameters
%%  lambda = 1e-1;
%%  lipschitz = 3;
%%  nIter = 5;
%%  
%%  % Initialise forward and inverse data
%%  b = sensor_data;
%%  x0 = p0_recon_adjoint.p_final;
%%  
%%  % Define operators
%%  forward_operator = @(initial_pressure) forward_operator_kWave(kgrid, medium, initial_pressure, sensor, input_args);
%%  inverse_operator = @(forward_data) inverse_operator_kWave(kgrid, medium, forward_data, sensor, input_args);
%%  %pixelKWave = fista(b, x0, lambda, lipschitz, nIter, forward_operator, inverse_operator);
%%  pixelKWave = gradient_descent(b, x0, lambda, nIter, forward_operator, inverse_operator);

%==================================================================================
% VISUALISATION
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex62_recon2D_subsample;

% Normalization
pixelRT = pixelRT./max(pixelRT(:));

% Axis
axisGrid = [0 1e3*Rgrid.xAxis(end) 0 1e3*Rgrid.yAxis(end)];

position = [700 700 320 600];
positionNoY = [700 700 300 600];
positionNoYBar = [700 700 363 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Initial Pressure
%==============================
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, Rgrid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gcf, 'pos', positionYBar);
saveas(gcf, 'Example62_U', 'png');
saveas(gcf, 'Example62_U.fig');
%title('Initial Pressure');

%========================================
% Initial reconstruction - RT
%========================================
pixelAReverse = Rgrid.pixelAReverse;
maxPixelRT = max(real(pixelAReverse(:)));
pixelRT_ini = max(0, real(pixelAReverse)/maxPixelRT);
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT_ini', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
set(gcf, 'pos', positionYBar);
%title('Reconstruction RT');
saveas(gcf, 'Example62_RT_subsample', 'png');
saveas(gcf, 'Example62_RT_subsample.fig');

%==============================
% Reconstruction - RT (FISTA)
%==============================
pixelRT = max(0, pixelRT_FISTA)/max(pixelRT_FISTA(:));
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x [mm]');
ylabel('y [mm]');
%set(gca, 'YTick', []);
set(gcf, 'pos', positionYBar);
%title('Reconstruction RT');
saveas(gcf, 'Example62_RT_FISTA_50iter', 'png');
saveas(gcf, 'Example62_RT_FISTA_50iter.fig');

%%  %==============================
%%  % Reconstruction - RT (SPDHG)
%%  %==============================
%%  pixelRT = max(0, pixelRT_SPDHG)/max(pixelRT_SPDHG(:));
%%  figure;
%%  surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, pixelRT', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  colorbar();
%%  box on;
%%  xlabel('x [mm]');
%%  ylabel('y [mm]');
%%  %set(gca, 'YTick', []);
%%  set(gcf, 'pos', positionYBar);
%%  %title('Reconstruction RT');
%%  saveas(gcf, 'Example62_RT_SPDHG_10epochs', 'png');
%%  saveas(gcf, 'Example62_RT_SPDHG_10epochs.fig');

%%  %========================================
%%  % Reconstruction - kWave
%%  %========================================
%%  pixelKWave = max(0, pixelKWave)/max(pixelKWave(:));
%%  figure;
%%  surf(grid.xAxis, grid.yAxis, pixelKWave', 'EdgeColor', 'none');
%%  view(2);
%%  axis(axisGrid);
%%  colorbar();
%%  box on;
%%  xlabel('x (m)');
%%  %ylabel('y (m)');
%%  set(gca, 'YTick', []);
%%  set(gcf, 'pos', positionNoYBar);
%%  %title('Reconstruction k-Wave (TR)');
%%  saveas(gcf, 'Example62_kWave_grad', 'png');

%========================================
% Error - RT
%========================================
%errorRT_TR = pixelRT - pixelKWave;
errorRT_TR = Rgrid.u - pixelRT;
figure;
surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, errorRT_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
%ylabel('y (m)');
set(gca, 'YTick', []);
set(gcf, 'pos', positionNoYBar);
%title('Error grad RT');
saveas(gcf, 'Example62_RT_errorGrad', 'png'); 

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);


