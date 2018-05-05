% Homogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex39_adjoint;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex38_gaussBeams;

%clear all;

%=========================================================================
% SIMULATION
%=========================================================================

%%  % create the computational grid
%%  Nx = 128;           % number of grid points in the x (row) direction
%%  Ny = 256;           % number of grid points in the y (column) direction
%%  dx = 1e-3;        % grid point spacing in the x direction [m]
%%  dy = 1e-3;        % grid point spacing in the y direction [m]
%%  kgrid = makeGrid(Nx, dx, Ny, dy);
%%  gridAux = gridRT(Nx, dx, Ny, dy);
%%  % define the properties of the propagation medium
%%  medium.sound_speed = 1500;
%%  medium.density = 1;
%%  
%%  % compute time
%%  dt = 3e-7;
%%  tMax = 4e-4;
%%  kgrid.t_array = 0:dt:tMax;
%%  %[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
%%  
%%  % Build initial pressure
%%  u0 = zeros(Nx, Ny);
%%  u0 = makeCircle(u0, 30, 128, 10, 1);
%%  u0 = makeCircle(u0, 60, 128, 10, 1);
%%  source.p0 = u0;
%%  % smooth the initial pressure distribution and restore the magnitude
%%  source.p0 = smooth(kgrid, source.p0, true);
%%  
%%  % input args
%%  input_args = {'PMLInside', false, 'Smooth', false};
%%  
%%  %==============================
%%  % Forward propagation
%%  %==============================
%%  
%%  % Define the sensors
%%  nSensors = 256;
%%  sensor.mask = zeros(Nx, Ny);
%%  sensor.mask(1, :) = 1;
%%  % run the simulation
%%  sensor.record = {'p', 'p_final'};
%%  sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  
%%  % Parallel sensors
%%  sensor.mask(:) = 0;
%%  sensor.mask(2, :) = 1;
%%  sensor_data_parallel = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  
%%  save sensor_data_forward.mat sensor sensor_data sensor_data_parallel kgrid source medium input_args dt;


%==============================
% Single sensor measurements
%==============================
%%  clear all;
    load sensor_data_forward.mat;
%%  
%%  sensor.mask(:) = 0;
%%  sensor.mask(10, 128) = 1;
%%  % Additive data
%%  source_single_add.p_mask = zeros(kgrid.Nx, kgrid.Ny);
%%  source_single_add.p_mask(1, :) = 1;
%%  source_single_add.p = fliplr(sensor_data.p);
%%  source_single_add.p_mode = 'additive';
%%  p0_recon_single_add = kspaceFirstOrder2D(kgrid, medium, source_single_add, sensor, input_args{:});
%%  
%%  % Dirichlet data
%%  source_single_dir.p_mask = zeros(kgrid.Nx, kgrid.Ny);
%%  source_single_dir.p_mask(1, :) = 1;
%%  source_single_dir.p = fliplr(sensor_data.p);
%%  source_single_dir.p_mode = 'dirichlet';
%%  p0_recon_dir = kspaceFirstOrder2D(kgrid, medium, source_single_dir, sensor, input_args{:});
%%  save recon_data_single_sensor.mat p0_recon_add p0_recon_dir;

%==============================
% Inverse propagation
%==============================
%%  %%%% Dirichlet derive data
%%  sensor_derive.mask = ones(kgrid.Nx, kgrid.Ny);
%%  sensor_derive.record = {'p_final'};
%%  sensor_data_derive = 0*sensor_data.p;
%%  sensor_data_derive(:, 2:end) = sensor_data.p(:, 2:end) - sensor_data.p(:, 1:end-1);
%%  source_derivative.p_mask = zeros(kgrid.Nx, kgrid.Ny);
%%  for n = 1:256
%%      disp(n);
%%      source_derivative.p_mask(:) = 0;
%%      source_derivative.p_mask(1, n) = 1;
%%      source_derivative.p = fliplr(sensor_data_derive(n, :));
%%      source_derivative.p_mode = 'dirichlet'; 
%%      p0_recon_deriv_dirichlet{n} = kspaceFirstOrder2D(kgrid, medium, source_derivative, sensor_derive, input_args{:});
%%  end
%%  save recon_data_dirichlet.mat p0_recon_deriv_dirichlet;

%%%% Dirichlet derive data
%%  sensor_der_normal.mask = ones(kgrid.Nx, kgrid.Ny);
%%  sensor_der_normal.record = {'p_final'};
%%  % Time derivative
%%  sensor_data_timeD = 0*sensor_data.p;
%%  sensor_data_timeD(:, 2:end-1) = (sensor_data.p(:, 3:end) - 2*sensor_data.p(:, 2:end-1) + sensor_data.p(:, 1:end-2))/dt/dt;
%%  sensor_data_timeD(:, 1) = (sensor_data.p(:, 3) - 2*sensor_data.p(:, 2) + sensor_data.p(:, 1))/dt/dt;
%%  sensor_data_timeD(:, end) = (sensor_data.p(:, end) - 2*sensor_data.p(:, end-1) + sensor_data.p(:, end-2))/dt/dt;
%%  % Space derivative
%%  sensor_data_spaceD = 0*sensor_data.p;
%%  sensor_data_spaceD(2:end-1, :) = (sensor_data.p(3:end, :) -2*sensor_data.p(2:end-1, :) + sensor_data.p(1:end-2, :))/(kgrid.dx/medium.sound_speed)^2;
%%  sensor_data_spaceD(1, :) = (sensor_data.p(3, :) -2*sensor_data.p(2, :) + sensor_data.p(1, :))/(kgrid.dx/medium.sound_speed)^2;
%%  sensor_data_spaceD(end, :) = (sensor_data.p(end, :) -2*sensor_data.p(end-1, :) + sensor_data.p(end-2, :))/(kgrid.dx/medium.sound_speed)^2;
% Sensor data
%%  sensor_data_der_normal = sensor_data.p - sensor_data_parallel.p;
%%  source_derivative.p_mask = zeros(kgrid.Nx, kgrid.Ny);
%%  for n = 1:256
%%      disp(n);
%%      source_derivative.p_mask(:) = 0;
%%      source_derivative.p_mask(1, n) = 1;
%%      source_derivative.p = fliplr(sensor_data_der_normal(n, :));
%%      source_derivative.p_mode = 'dirichlet'; 
%%      p0_recon_deriv_dirichlet_normal{n} = kspaceFirstOrder2D(kgrid, medium, source_derivative, sensor_der_normal, input_args{:});
%%  end
%%  save recon_data_derive_normal.mat p0_recon_deriv_dirichlet_normal;

%%  %%%% Additive data
%%  sensor.mask = ones(kgrid.Nx, kgrid.Ny);
%%  sensor.record = {'p_final'};
%%  source_additive.p_mask = zeros(kgrid.Nx, kgrid.Ny);
%%  for n = 1:256
%%      disp(n);
%%      source_additive.p_mask(:) = 0;
%%      source_additive.p_mask(1, n) = 1;
%%      source_additive.p = fliplr(sensor_data.p(n, :));
%%      source_additive.p_mode = 'additive'; 
%%      p0_recon_additive{n} = kspaceFirstOrder2D(kgrid, medium, source_additive, sensor, input_args{:});
%%  end
%%  save recon_data_additive.mat p0_recon_additive;

%%  %%%% Time reversal
%%  clear source;
%%  sensor.mask = ones(kgrid.Nx, kgrid.Ny);
%%  sensor.record = {'p_final'};
%%  source.p_mask = zeros(kgrid.Nx, kgrid.Ny);
%%  source.p_mask(1, :) = 1;
%%  source.p = fliplr(sensor_data.p);
%%  source.p_mode = 'dirichlet'; 
%%  p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%  save recon_data.mat p0_recon;


%=========================================================================
% VISUALISATION
%=========================================================================
close all;
%%  load sensor_data_forward.mat;
%%  load recon_data_single_sensor.mat;
%%  load recon_data_dirichlet.mat;
    load recon_data_additive.mat;
%%  load recon_data.mat;
%%  load recon_data_derive_normal.mat
%%  cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex39_adjoint;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex39_adjoint;

%%  %==============================
%%  % Sensor data
%%  %==============================
%%  colourMapV = cool(256);
%%  figure;
%%  hold on;
%%  for n = 1:256
%%      plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n, :));
%%  end
%%  ylabel('Sensor');
%%  xlabel('Time Step');
%%  saveas(gcf, 'Example39_kWave_signal.fig');

%%  %==============================
%%  % Signal for additive data
%%  %==============================
%%  figure;
%%  plot(kgrid.t_array, p0_recon_additive.p, 'Color', 'r');
%%  hold on;
%%  plot(kgrid.t_array, p0_recon_dirichlet.p, 'Color', 'b');

%%  %==============================
%%  % Signal for derivatives
%%  %==============================
%%  figure;
%%  plot(kgrid.t_array, sensor_data_timeD(1, :), 'Color', 'r');
%%  hold on;
%%  plot(kgrid.t_array, sensor_data_spaceD(1, :), 'Color', 'b');
%%  plot(kgrid.t_array, sensor_data_der_normal(1, :), 'Color', 'g');
%%  legend('time derivative', 'space derivative', 'difference');

%%  %==============================
%%  % Reconstruction for dirichlet
%%  %==============================
%%  p0_recon_dirichlet_sum = zeros(kgrid.Nx, kgrid.Ny);
%%  for n = 1:256
%%      p0_recon_dirichlet_sum = p0_recon_dirichlet_sum + p0_recon_deriv_dirichlet{n}.p_final;
%%  end
%%  %p0_recon_dirichlet_sum = max(0, p0_recon_dirichlet_sum);
%%  figure;
%%  surf(kgrid.x_vec, kgrid.y_vec, p0_recon_dirichlet_sum', 'EdgeColor', 'none');
%%  view(2);
%%  title('Dirichlet derivative')

%==============================
% Reconstruction for derivative
%==============================
p0_recon_derive_sum = zeros(kgrid.Nx, kgrid.Ny);
for n = 1:256
    p0_recon_derive_sum = p0_recon_derive_sum + p0_recon_deriv_dirichlet_normal{n}.p_final;
end
%p0_recon_dirichlet_sum = max(0, p0_recon_derive_sum);
figure;
surf(kgrid.x_vec, kgrid.y_vec, p0_recon_derive_sum', 'EdgeColor', 'none');
view(2);
title('Derivative using space')

%==============================
% Reconstruction for additive
%==============================
p0_recon_additive_sum = zeros(kgrid.Nx, kgrid.Ny);
for n = 1:256
    p0_recon_additive_sum = p0_recon_additive_sum + p0_recon_additive{n}.p_final;
end
%p0_recon_additive_sum = max(0, p0_recon_additive_sum);
figure;
surf(kgrid.x_vec, kgrid.y_vec, p0_recon_additive_sum', 'EdgeColor', 'none');
view(2);
title('Additive sources')

%==============================
% Reconstruction for difference
%==============================
p0_recon_both_sum = p0_recon_additive_sum - p0_recon_derive_sum;
p0_recon_both_sum = max(0, p0_recon_both_sum);
figure;
surf(kgrid.x_vec, kgrid.y_vec, p0_recon_both_sum', 'EdgeColor', 'none');
view(2);
title('Sum additive dirichlet')

%%  %==============================
%%  % Reconstruction time reversal
%%  %==============================
%%  p0_recon.p_final = max(0, p0_recon.p_final);
%%  figure;
%%  surf(kgrid.x_vec, kgrid.y_vec, p0_recon.p_final', 'EdgeColor', 'none');
%%  view(2);
%%  title('Time reversal')

%%  %==============================
%%  % Error
%%  %==============================
%%  error = p0_recon_both_sum/max(p0_recon_both_sum(:)) - p0_recon.p_final/max(p0_recon.p_final(:));
%%  figure;
%%  surf(kgrid.x_vec, kgrid.y_vec, error', 'EdgeColor', 'none');
%%  view(2);
%%  title('Error')

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
