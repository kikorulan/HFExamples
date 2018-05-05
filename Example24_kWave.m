% Homogeneous Propagation Medium Example

clear all;

run colourMap;
% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 200;           % number of grid points in the x (row) direction
Ny = 200;           % number of grid points in the y (column) direction (multiple of 3)
dx = 0.1/Nx;        % grid point spacing in the x direction [m]
dy = 0.1/Ny;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);

%%%%%%%%%%%%%%%%%%
% define the properties of the propagation medium
%%%%%%%%%%%%%%%%%
% Build domain
c0 = 1500;
% Randomise C
c = c0*ones(Nx, Ny);
gridAux.setCMatrix(c);
spatialSDV = 0.09;
deltaC = 0.1;
gridAux.randomiseC(spatialSDV, deltaC);

medium.sound_speed = gridAux.c;
medium.density = 1;

% compute time
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/4),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M2 ; M3] M4];

source.p0 = U;

% Define the sensors
nSensors = 3;
% Source 1
x1 = [Nx*dx/2; Ny*dy/100]; 
point = gridAux.findCoordinates(x1);
C1 = gridAux.getC(point);
gridAux.setU(point, 1);
% Source 2
x2 = [5*Nx*dx/6; Ny*dy/100];
point = gridAux.findCoordinates(x2);
C2 = gridAux.getC(point);
gridAux.setU(point, 1);
% Source 3
x3 = [Nx*dx/100; Ny*dy/3]; 
point = gridAux.findCoordinates(x3);
C3 = gridAux.getC(point);
gridAux.setU(point, 1);

sensor.mask = gridAux.u;
size(gridAux.u);
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
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_RT_randomC/Example24_U_kWave.fig');


% plot the simulated sensor data
figure;
hold on;
for n = 1:nSensors
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n));
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex24_RT_randomC/Example24_aSignal_kWave.fig');
save('workspace');
