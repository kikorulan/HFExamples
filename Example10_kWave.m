% Homogeneous Propagation Medium Example

clear kgrid source;

run colorMap;
% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 1080;           % number of grid points in the x (row) direction
Ny = 1080;           % number of grid points in the y (column) direction
dx = 0.12/Nx;        % grid point spacing in the x direction [m]
dy = 0.12/Ny;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);


% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
% compute time
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);


% create initial pressure distribution using makeDisc
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];

source.p0 = U;

% Define the sensors
nSensors = 3;
x1 = [Nx*dx/2; Ny*dy/100]; % Source 1
x2 = [5*Nx*dx/6; Ny*dy/100]; % Source 2
x3 = [Nx*dx/100; Ny*dy/3]; % Source 3
point = gridAux.findCoordinates(x1);
gridAux.setU(point, 1);
point = gridAux.findCoordinates(x2);
gridAux.setU(point, 1);
point = gridAux.findCoordinates(x3);
gridAux.setU(point, 1);
sensor.mask = gridAux.u;

% run the simulation
sensor.record = {'p', 'p_final'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, fliplr(source.p0)' + fliplr(sensor.mask)', [-1 1]);
colormap(getColorMap);
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_U_kWave.fig');


% plot the simulated sensor data
figure;
hold on;
for n = 1:nSensors
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMap(n));
    plot(t1, 10*s1(n, :), 'Color', colourMap(n+nSensors));
end
legend('Sensor 1', 'Sensor 1 - 1e-3', 'Sensor 2', 'Sensor 2 - 1e-3', 'Sensor 3', 'Sensor 3 - 1e-3');
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_aSignal_kWave.fig');
