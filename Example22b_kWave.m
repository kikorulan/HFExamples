% Homogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;

clear all;

run colourMap;
% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 240;%            % number of grid points in the x (row) direction
Ny = 360;%           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
gridAux = gridRT(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500;
v1 = -0.1;
v2 = 0.1;
kernelSize = 33;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, floor(dimY/3));
M2 = v1*c0*ones(floor(dimX/2), floor(dimY/3));
M3 = v2*c0*ones(floor(dimX/2), floor(dimY/3));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [M1...
     [M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
cMatrix = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.sound_speed = cMatrix;
medium.density = 1;

% compute time
dt = 3e-8;
tMax = 4e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);


% create initial pressure distribution
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];

source.p0 = U;
source.p0 = smooth(kgrid, source.p0, true);

% Define the sensors
nSensors = 3;
% Sources locations
x1(1, 1, 1) = Nx*dx/2; 
x1(1, 1, 2) = Ny*dy/100; 
x2(1, 1, 1) = 5*Nx*dx/6; 
x2(1, 1, 2) = Ny*dy/100; 
x3(1, 1, 1) = Nx*dx/100; 
x3(1, 1, 2) = Ny*dy/2; 
point = gridAux.findCoordinates(x1);
gridAux.setU(point, 1);
point = gridAux.findCoordinates(x2);
gridAux.setU(point, 1);
point = gridAux.findCoordinates(x3);
gridAux.setU(point, 1);
sensor.mask = gridAux.u;
size(gridAux.u)
% run the simulation
sensor.record = {'p', 'p_final'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);

save kgrid_data.mat sensor_data kgrid source cMatrix c0;

% =========================================================================
% VISUALISATION
% =========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex22_kWave_interf;


% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, fliplr(source.p0)' + fliplr(sensor.mask)', [-1 1]);
colormap(getColorMap);
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
saveas(gcf, 'Example22b_U_kWave.fig');

% plot the simulated sensor data
figure;
hold on;
for n = 1:nSensors
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n));
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, 'Example22b_aSignal_kWave.fig');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/;
