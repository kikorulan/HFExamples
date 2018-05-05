% Homogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

clear all;
close all;

run colourMap;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
medium.sound_speed = c0;
medium.density = 1;

% compute time
dt = 3e-7;
tMax = 2e-4;
kgrid.t_array = 0:dt:tMax;

% Build initial pressure
u0 = zeros(Nx, Ny);
u0 = makeCircle(u0, 30, 128, 10, 1);
u0 = makeCircle(u0, 60, 128, 10, 1);
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);

%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
%sensor.mask(1, :) = 1;
sensor.mask(1, :) = 1;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

save sensor_data.mat sensor_data kgrid source sensor medium u0;
% reset the initial pressure
source.p0 = 0;

%==============================
% Time reversal
%==============================
%%  % assign the time reversal data
%%  sensor.time_reversal_boundary_data = sensor_data;
%%  clear p0_recon;
%%  for n = 1:256
%%      disp(n);
%%      sensor.mask(:) = 0;
%%      sensor.mask(1, n) = 1;
%%      sensor.time_reversal_boundary_data = sensor_data(n, :);
%%      p0_recon{n} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  end
%%  sensor.mask(1, :) = 1;
%%  sensor.time_reversal_boundary_data = sensor_data;
%%  p0_recon{257} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  save recon_data_sensors.mat p0_recon kgrid;

%==============================
% Time varying source
%==============================
clear p0_recon_TVsource source;
sensor.mask(:) = 1;
sensor.record = {'p_final'};
source.p_mask = zeros(Nx, Ny);
for n = 1:256
    disp(n);
    source.p_mask(:) = 0;
    source.p_mask(1, n) = 1;
    source.p = fliplr(sensor_data(n, :));
    p0_recon_TVsource{n} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end
source.p_mask(1, :) = 1;
source.p = fliplr(sensor_data);
p0_recon_TVsource{257} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
save recon_data_sensors_TVsource.mat p0_recon_TVsource kgrid;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

% plot the initial pressure and sensor distribution
%%  figure;
%%  imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, fliplr(source.p0)' + fliplr(sensor.mask)', [-1 1]);
%%  colormap(getColorMap);
%%  ylabel('y-position [mm]');
%%  xlabel('x-position [mm]');
%%  axis image;
%%  saveas(gcf, 'Example31_U_kWave.fig');

% plot the simulated sensor data
%%  figure;
%%  hold on;
%%  for n = 1:nSources
%%      plot(kgrid.t_array, sensor_data(n, :), 'Color', colourMapV(n));
%%  end
%%  legend('Sensor 1', 'Sensor 2', 'Sensor 3');
%%  ylabel('Sensor');
%%  xlabel('Time Step');
%%  saveas(gcf, 'Example31_aSignal_kWave.fig');

% plot the reconstructed initial pressure 
%%  figure;
%%  imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, p_xy_rs, [-1 1]);
%%  colormap(getColorMap);
%%  ylabel('x-position [mm]');
%%  xlabel('y-position [mm]');
%%  axis image;
%%  colorbar();
%%  saveas(gcf, 'Example31_Urecon_kWave');


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
