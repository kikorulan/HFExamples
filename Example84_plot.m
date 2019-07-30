% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex84_ROI
clear all;
close all;


% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%==================================================
% Dimensions
%==================================================
% Import dimensions
filenameDimensions = './input_data/dimensions.dat';
dim = importdata(filenameDimensions, delimiterIn, headerlinesIn);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%================================================================================
% INITIAL PRESSURE AND ROI
%================================================================================
% Initial Pressure
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
a = axes;
t = title('Initial Pressure');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_initialPressure.fig')

% Scatter ROI
sensor_roi = importdata('./input_data/sensors_ROI_3600.dat', delimiterIn, headerlinesIn);
xpos = sensor_roi(1, 1);
ypos = sensor_roi(1, 2);
zpos = sensor_roi(1, 3);
radius = sensor_roi(1, 4);
scatter_points = scatter_roi(Nx, Ny, Nz, dx, dy, dz, radius, xpos, ypos, zpos);
f = figure;
plot3(scatter_points(:, 2), scatter_points(:, 3), scatter_points(:, 1), 'or');
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Ny*dy 0 Nz*dz 0 Nx*dx]);
daspect([1 1 1]);
f.CurrentAxes.ZDir = 'Reverse';
title('Region of Interest');
saveas(gcf, './figures/Example84_ROI.fig')

% Initial Pressure ROI
u_roi = select_roi(u0, dx, dy, dz, radius, xpos, ypos, zpos);
h = plot_projection(u_roi, dx);
a = axes;
t = title('Initial Pressure (ROI mask)');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_initialPressure_ROI.fig')

%================================================================================
% PLOT ROI
%================================================================================
% Single BEAM
single_trajectory = importdata(['./output_data/Trajectory5.dat'], delimiterIn, headerlinesIn);
% Read number of rays and steps
[nSteps single_nRays] = size(single_trajectory);
single_xCoord = single_trajectory(:, 1:3:single_nRays);
single_yCoord = single_trajectory(:, 2:3:single_nRays);
single_zCoord = single_trajectory(:, 3:3:single_nRays);
single_nRays = single_nRays/3;

% Plot
f = figure;
plot3(scatter_points(:, 2), scatter_points(:, 3), scatter_points(:, 1), 'or');
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Ny*dy 0 Nz*dz 0 Nx*dx]);
daspect([1 1 1]);
f.CurrentAxes.ZDir = 'Reverse';
hold on;
colours = winter(single_nRays);
for n = 1:5:single_nRays
    plot3(single_yCoord(:, n), single_zCoord(:, n), single_xCoord(:, n), 'Color', colours(n, :));
end
title('ROI illuminated by single beam');
saveas(gcf, './figures/Example84_beam', 'png')

% Random Trajectories
nTrajectories = 5;
nSensors = 36;
index_trajec = 1:nSensors;
rng(1);
index_random = index_trajec(randperm(length(index_trajec)));
for i = 1:nTrajectories
    index = index_random(i);
    % Import data
    trajectories{i} = importdata(['./output_data/Trajectory' int2str(index-1) '.dat'], delimiterIn, headerlinesIn);
    % Read number of rays and steps
    [nSteps nRays{i}] = size(trajectories{i});
    xCoord{i} = trajectories{i}(:, 1:3:nRays{i});
    yCoord{i} = trajectories{i}(:, 2:3:nRays{i});
    zCoord{i} = trajectories{i}(:, 3:3:nRays{i});
    nRays{i} = nRays{i}/3;
end

% Plot
f = figure;
plot3(scatter_points(:, 2), scatter_points(:, 3), scatter_points(:, 1), 'or');
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Ny*dy 0 Nz*dz 0 Nx*dx]);
daspect([1 1 1]);
f.CurrentAxes.ZDir = 'Reverse';
colours = winter(nTrajectories);
hold on;
for i = 1:nTrajectories
    for n = 1:10:nRays{i}
        plot3(yCoord{i}(:, n), zCoord{i}(:, n), xCoord{i}(:, n), 'Color', colours(i, :));
    end
end
title('ROI illuminated by multiple beams');
saveas(gcf, './figures/Example84_multiBeam', 'png')

%================================================================================
% PLOT FULL ADJOINT
%================================================================================
pressureFullAdj_Matrix = importdata('./input_data/pixelPressure_adjoint.dat', ' ', 0);
pressureFullAdj = matrix2cube(pressureFullAdj_Matrix, Nz);
h = plot_projection(pressureFullAdj, dx);
a = axes;
t = title('Full Adjoint Pressure');
a.Visible = 'off'; 
t.Visible = 'on'; 

pressureFullAdj_roi = select_roi(pressureFullAdj, dx, dy, dz, radius, xpos, ypos, zpos);
h = plot_projection(pressureFullAdj_roi, dx);
a = axes;
t = title('Full Adjoint Pressure ROI');
a.Visible = 'off'; 
t.Visible = 'on'; 

%================================================================================
% PLOT ROI ADJOINT
%================================================================================
pressureAdj_Matrix = importdata('./output_data/pixelPressure_adjoint_ROI.dat', ' ', 0);
pressureAdj = matrix2cube(pressureAdj_Matrix, Nz);
h = plot_projection(pressureAdj, dx);
a = axes;
t = title('Adjoint Pressure ROI operator');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_adjoint.fig')

pressureAdj_roi = select_roi(pressureAdj, dx, dy, dz, radius, xpos, ypos, zpos);
h = plot_projection(pressureAdj_roi, dx);
a = axes;
t = title('Adjoint Pressure ROI operator (ROI mask)');
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_adjoint_ROI.fig')

%==================================================================================================================================
% PLOT ADJOINT
%==================================================================================================================================
%==============================
% Gradient Descent
%==============================
GD = [];
GD.tau    = '1e18';
GD.lambda = '1e-4';
GD.iter   = int2str(10);
pixelPressureMatrix = importdata(['./output_data/pixelPressure_GD_tau', GD.tau, '_lambda', GD.lambda, '_iter', GD.iter, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['Gradient Descent - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter]);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_reconGD_10iter.fig')

pixelPressure_roi = select_roi(pixelPressure, dx, dy, dz, radius, xpos, ypos, zpos);
h = plot_projection(pixelPressure_roi, dx);
a = axes;
t = title(['Gradient Descent - t = ', GD.tau, ', l = ', GD.lambda, ', iter = ', GD.iter, ' (ROI mask)']);
a.Visible = 'off'; 
t.Visible = 'on';
saveas(gcf, './figures/Example84_reconGD_10iter_ROI.fig')
 
%==============================
% Stochastic Gradient Descent
%==============================
SGD = [];
SGD.tau    = '4e18';
SGD.lambda = '1e-4';
SGD.batch  = '100';
SGD.epoch  = int2str(10);
pixelPressureMatrix = importdata(['./output_data/pixelPressure_S-GD_tau', SGD.tau, '_lambda', SGD.lambda, '_batch', SGD.batch, '_subepoch', SGD.epoch, '.dat'], ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);
a = axes;
t = title(['Stochastic GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', epoch = ', SGD.epoch]);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_reconSGD_10epoch.fig')

pixelPressure_roi = select_roi(pixelPressure, dx, dy, dz, radius, xpos, ypos, zpos);
h = plot_projection(pixelPressure_roi, dx);
a = axes;
t = title(['Stochastic GD - t = ', SGD.tau, ', l = ', SGD.lambda, ', epoch = ', SGD.epoch, ' (ROI mask)']);
a.Visible = 'off'; 
t.Visible = 'on'; 
saveas(gcf, './figures/Example84_reconSGD_10epoch_ROI.fig')
