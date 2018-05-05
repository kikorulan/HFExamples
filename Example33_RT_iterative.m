%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex33_reconstruction;

close all;
% Measure computational time
tic;
start_time = clock;

%load sensor_data.mat;
%load recon_data_sensors.mat;
load recon_data_adjoint.mat;
load recon_data_RT.mat;
%==================================================================================
% SIMULATION
%==================================================================================
%==============================
% RT
%==============================
% Define parameters
%%  lambda = 1e-1;
%%  lipschitz = 3;
%%  nIter = 5;
lambda = 1;
sigma = 1e4;
tau = 1e8;
theta = 1;
nEpochs = 1;
nSensors = 256;

% Parameters
para.maxIter = 20;
para.constraint = 'positivity';

% Define operators
%norm_factor = max(grid.u(:))/max(grid.pixelAReverse(:));
norm_factor = 1;
forward_operator = @(index, initial_pressure) grid.operator_forward(source, index, initial_pressure, norm_factor);
inverse_operator = @(index, forward_data) grid.operator_inverse(source, index, forward_data);

% Initialise forward and inverse data
if (exist('y0')); clear y0; end;
for n = 1:length(source)
    if(isempty(source(n).aForward_initial))
        source(n).setForwardSignal_initial(source(n).aForward);
    end
    y0(n, :) = source(n).aForward_initial;
end
x0 = pixelAReverse;

% Proximal operators
proxG = @(x, a) conTVdenoising(x, 1e-2, para);
proxF = @(y, a, b) (y-a*y0(b, :))/(1+a); 

% Run iterations
%pixelRT = fista(b, x0, lambda, lipschitz, nIter, forward_operator, inverse_operator);
%pixelRT = gradient_descent(b, norm_factor*x0, lambda, nIter, forward_operator, inverse_operator);
pixelRT = StochasticPDHG(x0, y0, proxF, proxG, forward_operator, inverse_operator, sigma, tau, theta, nEpochs, nSensors);
pixelRT = pixelRT./max(pixelRT(:));
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
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex33_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex33_reconstruction;
axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 320 600];
positionNoY = [700 700 300 600];
positionNoYBar = [700 700 363 600];
positionYBar = [700 700 390 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Initial Pressure
%==============================
figure;
surf(grid.xAxis, grid.yAxis, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
%title('Initial Pressure');
saveas(gcf, 'Example33_U', 'png'); 

%==============================
% Reconstruction - RT
%==============================
pixelRT = max(0, pixelRT)/max(pixelRT(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
%colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
%set(gca, 'YTick', []);
set(gcf, 'pos', position);
%title('Reconstruction RT');
saveas(gcf, 'Example33_RT_grad', 'png');

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
%%  saveas(gcf, 'Example33_kWave_grad', 'png');

%========================================
% Error - RT
%========================================
%errorRT_TR = pixelRT - pixelKWave;
errorRT_TR = grid.u - pixelRT;
figure;
surf(grid.xAxis, grid.yAxis, errorRT_TR', 'EdgeColor', 'none');
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
saveas(gcf, 'Example33_RT_errorGrad', 'png'); 

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

