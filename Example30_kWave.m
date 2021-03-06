% Homogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex30_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex30_interf;

clear kgrid source;

run colourMap;
% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
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
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.density = 1;

% compute time
dt = 3e-7;
tMax = 4e-4;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);


% create initial pressure distribution
M0 = ones(floor(Nx/4), floor(Ny/6));
M1 = 0*M0;
M2 = triu(M0);
M3 = fliplr(triu(M0)')';
U = [M1 M1 M1 M1 M1 M1;
    M1 M3 M2' M1 M1 M1;
    M1 M2 fliplr(M2) M1 M1 M1;
    M1 M1 M1 M1 M1 M1];

sourceKW.p0 = U;
% smooth the initial pressure distribution and restore the magnitude
sourceKW.p0 = smooth(kgrid, sourceKW.p0, true);

% Define the sensors
nSensors = 3;
x1(1, 1, 1) = 2*Nx*dx/3;
x1(1, 1, 2) = 0; % Source 1
x2(1, 1, 1) = Nx*dx/4;
x2(1, 1, 2) = (Ny-1)*dy; % Source 2
x3(1, 1, 1) = 3*Nx*dx/4; 
x3(1, 1, 2) = (Ny-1)*dy; % Source 3
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
sensor_data = kspaceFirstOrder2D(kgrid, medium, sourceKW, sensor, 'PMLInside', false);

save sensor_data.mat sensor_data sourceKW kgrid medium;

% =========================================================================
% VISUALISATION
% =========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex30_interf;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex30_interf;

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, fliplr(sourceKW.p0)' + fliplr(sensor.mask)', [-1 1]);
colormap(getColorMap);
ylabel('y-position [mm]');
xlabel('x-position [mm]');
axis image;
saveas(gcf, 'Example30_U_kWave.fig');

% plot the simulated sensor data
figure;
hold on;
for n = 1:nSensors
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n));
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, 'Example30_aSignal_kWave.fig');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
