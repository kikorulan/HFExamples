% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex75_3D_thinveins_het;

clear all;
close all;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================
x_axis = 1e3*(0:dx:(Nx-1)*dx);
y_axis = 1e3*(0:dy:(Ny-1)*dy);
z_axis = 1e3*(0:dz:(Nz-1)*dz);
ssLim = 1580*[0.95 1.05];
% Position
position     = [700 700 300 630];
positionY    = [700 700 320 630];
positionBar  = [700 700 363 630];
positionYBar = [700 700 390 630];
set(0,'DefaultFigurePaperPositionMode','auto');

%==================================================
% SOUND SPEED
%==================================================
ssFigPosition = [700 700 1500 300];
ssAxisPosition1 = [0.05 0.1 0.40 0.9];
ssAxisPosition2 = [0.5 0.1 0.40 0.9];
cMatrix = importdata('input_data/sound_speed.dat', ' ', 0);
c = matrix2cube(cMatrix, Nz);

figure;
subplot(1, 2, 1);
imagesc(y_axis, x_axis, c(:, :, 80)');
xlabel('y [mm]');
ylabel('x [mm]');
caxis(ssLim);
pbaspect([3 1 1]);
set(gca, 'position', ssAxisPosition1);

subplot(1, 2, 2);
imagesc(y_axis, x_axis, c(:, :, 160)');
xlabel('y [mm]');
ylabel('x [mm]');
caxis(ssLim);
pbaspect([3 1 1]);
colorbar();
set(gca, 'position', ssAxisPosition2);
set(gcf, 'position', ssFigPosition);
saveas(gcf, './figures/pID_sound_speed.fig');
saveas(gcf, './figures/pID_sound_speed', 'png');
%==================================================
% INITIAL PRESSURE
%==================================================
% Load Initial Pressure
u0Matrix = importdata('input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
saveas(gcf, './figures/pID_initial_pressure.fig');
saveas(gcf, './figures/pID_initial_pressure', 'png');
% Load Initial Pressure
load ./input_data/initial_pressure_veins_smooth;
plot_projection(initial_pressure_veins_smooth, dx);

%==================================================
% Reconstruction - kWave
%==================================================
PML_size = 10;
% 14400 sensors
p0_recon_PML = h5read('output_data/Example75_adjoint_output_14400sensors.h5', '/p_final');
p0_recon = p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
plot_projection(p0_recon, dx);
saveas(gcf, './figures/pID_kWave_adjoint_14400sensors.fig');
saveas(gcf, './figures/pID_kWave_adjoint_14400sensors', 'png');

% 1600 sensors
p0_recon_PML = h5read('output_data/Example75_adjoint_output_1600sensors.h5', '/p_final');
p0_recon = p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
plot_projection(p0_recon, dx);
saveas(gcf, './figures/pID_kWave_adjoint_1600sensors.fig');
saveas(gcf, './figures/pID_kWave_adjoint_1600sensors', 'png');

%==================================================
% ADJOINT RT
%==================================================
% 14400
pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_RT_14400sensors.dat', ' ', 0);
pixelPressure = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure, dx);
pixelPressureMatrix_norm = pixelPressureMatrix*normRT;
dlmwrite('input_data/pixelPressure_adjoint_RT_14400sensors_norm1e18.dat', pixelPressureMatrix_norm, 'delimiter', ' ');

% 1600
pixelPressureMatrix = importdata('input_data/pixelPressure_adjoint_RT_1600sensors.dat', ' ', 0);
pixelPressure = matrix2cube(pixelPressureMatrix, Nz);
plot_projection(pixelPressure, dx);

%==================================================
% ITERATIVE RECONSTRUCTION
%==================================================
for iter = 1:5
    % Gradient
    pixelPressureMatrix = importdata(['results/pixelPressure_GD_tau2e18_lambda1e-2_iter', int2str(iter), '.dat'], ' ', 0);
    pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    plot_projection(pixelPressure, dx);
    a = axes;
    t = title(['GD - t = 1e18, l = 1e-2, iter = ', int2str(iter), ' - homogeneous SS']);
    a.Visible = 'off'; 
    t.Visible = 'on'; 
    %saveas(gcf, ['figures/Example75_GD_tau1e18_lambda1e-2_iter', int2str(iter), '.fig']);
    
    % Stochastic gradient
    pixelPressureMatrix = importdata(['results/pixelPressure_S-GD_tau8e18_lambda3e-4_batch90_epoch', int2str(iter), '.dat'], ' ', 0);
    pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    plot_projection(pixelPressure, dx);
    a = axes;
    t = title(['S-GD - t = 4e18, l = 3e-4, batch = 90, epoch = ', int2str(iter), ' - homogeneous SS']);
    a.Visible = 'off'; 
    t.Visible = 'on'; 
    %saveas(gcf, ['figures/Example75_S-GD_tau2e18_lambda3e-4_batch90_epoch', int2str(iter), '.fig']);
    
    % FISTA
    pixelPressureMatrix = importdata(['results/pixelPressure_FISTA_tau1e18_lambda1e-2_iter', int2str(iter), '.dat'], ' ', 0);
    pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    plot_projection(pixelPressure, dx);
    a = axes;
    t = title(['FISTA - t = 1e18, l = 1e-2, iter = ', int2str(iter), ' - homogeneous SS']);
    a.Visible = 'off'; 
    t.Visible = 'on'; 
    %saveas(gcf, ['figures/Example75_FISTA_tau1e18_lambda1e-2_iter', int2str(iter), '.fig']);
    
    %%  % Stochastic FISTA
    %%  pixelPressureMatrix = importdata('results/pixelPressure_S-FISTA_tau1e18_lambda5e-4_epoch50.dat', ' ', 0);
    %%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    %%  plot_projection(pixelPressure, dx);
    %%  a = axes;
    %%  t = title('GD - t = 1e18, l = 1e-2, iter = 50 - homogeneous SS');
    %%  a.Visible = 'off'; 
    %%  t.Visible = 'on'; 
    %%  saveas(gcf, 'Example68_GD_tau1e18_lambda1e-2_iter50.fig');
    
    %%  % PDHG
    %%  pixelPressureMatrix = importdata('results/pixelPressure_PDHG_sigma1e0_tau1e18_theta1_lambda1e-2_iter1.dat', ' ', 0);
    %%  pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    %%  plot_projection(pixelPressure, dx);
    %%  a = axes;
    %%  t = title('PDHG - s = 1, t = 1e18, l = 1e-2, iter = 50 - homogeneous SS');
    %%  a.Visible = 'off'; 
    %%  t.Visible = 'on'; 
    %%  %saveas(gcf, 'figures/Example75_PDHG_sigma1_tau1e18_lambda1e-2_iter50.fig');
    
    % S-PDHG
    pixelPressureMatrix = importdata(['results/pixelPressure_S-PDHG_sigma2e1_tau8e18_theta1_lambda3e-4_batch90_epoch', int2str(iter), '.dat'], ' ', 0);
    pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
    plot_projection(pixelPressure, dx);
    a = axes;
    t = title(['S-PDHG - s = 2e0, t = 4e18, l = 3e-4, epoch = ', int2str(iter), ' - homogeneous SS']);
    a.Visible = 'off'; 
    t.Visible = 'on'; 
    %saveas(gcf, ['figures/Example75_PDHG_sigma2e0_tau4e18_lambda1e-2_epoch', int2str(iter), '.fig']);
end

%==================================================
% AUXILIAR RECONSTRUCTION
%==================================================
% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_GD_2.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_S-GD_2.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_FISTA_0.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_S-FISTA_29.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_PDHG_1.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('./output_data/pixelPressure_S-PDHG_34.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);

% Import data
pixelPressureMatrix = importdata('output_data/pixelPressure_TVdenoised_1e-14-3000.dat', ' ', 0);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure, dx);




difference = u0 - pixelPressure;
plot_projection(difference, dx);
