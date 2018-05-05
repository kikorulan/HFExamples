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

% Single sensor 100
nSensor = 100;
source.p_mask(:) = 0;
source.p_mask(1, nSensor) = 1;
source.p = fliplr(sensor_data(nSensor, :));
source.p_mode = 'dirichlet';
p0_recon_TVsource{1} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
source.p_mode = 'additive';
p0_recon_TVsource{2} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Single sensor 200
nSensor = 200;
source.p_mask(:) = 0;
source.p_mask(1, nSensor) = 1;
source.p = fliplr(sensor_data(nSensor, :));
source.p_mode = 'dirichlet';
p0_recon_TVsource{3} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
source.p_mode = 'additive';
p0_recon_TVsource{4} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Two sensors: 100 & 200
nSensor = [100 200];
source.p_mask(:) = 0;
source.p_mask(1, nSensor) = 1;
source.p = fliplr(sensor_data(nSensor, :));
source.p_mode = 'dirichlet';
p0_recon_TVsource{5} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
source.p_mode = 'additive';
p0_recon_TVsource{6} = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

save recon_data_dirichlet.mat p0_recon_TVsource kgrid;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

X = gridAux.xAxis;
Y = gridAux.yAxis;
%==============================
% Time varying source
%==============================
% Single sensor 100 - Dirichlet
figure;
surf(X, Y, p0_recon_TVsource{1}.p_final', 'EdgeColor', 'none'); 
view(2);
title('Single sensor 100- Dirichlet');

% Single sensor 100 - Additive
figure;
surf(X, Y, p0_recon_TVsource{2}.p_final', 'EdgeColor', 'none'); 
view(2);
title('Single sensor 100 - Additive');

% Single sensor 200 - Dirichlet
figure;
surf(X, Y, p0_recon_TVsource{3}.p_final', 'EdgeColor', 'none'); 
view(2);
title('Single sensor 200- Dirichlet');

% Single sensor 200 - Additive
figure;
surf(X, Y, p0_recon_TVsource{4}.p_final', 'EdgeColor', 'none'); 
view(2);
title('Single sensor 200 - Additive');

% Two sensors Dirichlet
figure;
surf(X, Y, p0_recon_TVsource{5}.p_final', 'EdgeColor', 'none'); 
view(2);
title('Two sensors - Dirichlet');

% Two sensors Additive
figure;
surf(X, Y, p0_recon_TVsource{6}.p_final', 'EdgeColor', 'none'); 
view(2);
title('Two sensors - Additive');

%==============
% Errors
%==============
% Sum Dirichlet
sumDirichlet = p0_recon_TVsource{1}.p_final + p0_recon_TVsource{3}.p_final;
figure;
surf(X, Y, sumDirichlet', 'EdgeColor', 'none'); 
view(2);
colorbar();
title('Sum Dirichlet');

% Sum Additive
sumAdditive = p0_recon_TVsource{2}.p_final + p0_recon_TVsource{4}.p_final;
figure;
surf(X, Y, sumAdditive', 'EdgeColor', 'none'); 
view(2);
colorbar();
title('Sum Additive');

% Error Dirichlet
errorDirichlet = p0_recon_TVsource{5}.p_final - sumDirichlet;
figure;
surf(X, Y, errorDirichlet', 'EdgeColor', 'none'); 
view(2);
colorbar();
title('Error Dirichlet');

% Error Additive
errorAdditive = p0_recon_TVsource{6}.p_final - sumAdditive;
figure;
surf(X, Y, errorAdditive', 'EdgeColor', 'none'); 
view(2);
colorbar();
title('Error Additive');


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
