% Sensor Element Directivity in 2D Example

clear all
close all

% =========================================================================
% DEFINE THE GRID AND MEDIUM PROPERTIES
% =========================================================================

%%  % create the computational grid
%%  Nx = 128;        % number of grid points in the x (row) direction
%%  Ny = 128;        % number of grid points in the y (column) direction
%%  dx = 1e-3/Nx;   % grid point spacing in the x direction [m]
%%  dy = dx;     	% grid point spacing in the y direction [m]
%%  kgrid = makeGrid(Nx, dx, Ny, dy);
%%  
%%  % define the properties of the propagation medium
%%  medium.sound_speed = 1500;  % [m/s]
%%  
%%  % create the time array, then halve it to make the example end sooner
%%  [kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
%%  kgrid.t_array = (0:0.7*length(kgrid.t_array))*dt;
%%  
%%  % =========================================================================
%%  % DEFINE THE DIRECTIONAL SENSOR ARRAY
%%  % =========================================================================
%%  
%%  % define a line of sensor points
%%  sensor.mask = zeros(Nx,Ny);
%%  sensor.mask(32,2:2:127) = 1;
%%  
%%  % define the angle of max directivity for each sensor point:
%%  %    0             = max sensitivity in x direction (up/down)
%%  %    pi/2 or -pi/2 = max sensitivity in y direction (left/right)
%%  sensor.directivity_angle = zeros(Nx,Ny);
%%  sensor.directivity_angle(32,2:2:127) = (-2:4/62:2)*pi/2;
%%  
%%  % define the directivity pattern
%%  sensor.directivity_pattern = 'pressure';
%%  
%%  % define the directivity size
%%  sensor.directivity_size = 16*kgrid.dx;
%%  
%%  % =========================================================================
%%  % SIMULATION AND VISUALISATION FOR AN INITIAL VALUE PROBLEM
%%  % =========================================================================
%%  
%%  % define the initial pressure distribution
%%  source.p0 = zeros(Nx,Ny);
%%  source.p0(63:65, :) = 2;
%%   
%%  % run the simulation
%%  sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLAlpha', [2 0]);
%%  
%%  % plot the largest value of the output for each sensor
%%  figure;
%%  plot(((-2:4/62:2)*pi/2), max(sensor_data, [], 2), 'o');
%%  xlabel('sensor directivity angle (radians)')
%%  ylabel('maxima of single-element sensors'' outputs')
%%  
%%  figure;
%%  max_data = max(sensor_data(:, 150:end), [], 2);
%%  polar(((-2:4/62:2)'*pi/2), max(sensor_data, [], 2));
