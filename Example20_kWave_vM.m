% Homogeneous Propagation Medium Example

clear all;
close all;
% Define colours
colourMapV = containers.Map('KeyType', 'double', 'ValueType', 'any');
colourMapV(1) = [1 0 0];
colourMapV(2) = [0 1 0];
colourMapV(3) = [0 0 1];

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 120;           % number of grid points in the x (row) direction
Ny = 180;           % number of grid points in the y (column) direction
dx = 0.12/Nx;        % grid point spacing in the x direction [m]
dy = 0.18/Ny;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

%%%%%%%%%%%%%%%%%%
% define the properties of the propagation medium
%%%%%%%%%%%%%%%%%
% Build domain
c0 = 1500;
v1 = 5*-0.1;
v2 = 5*0.1;
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
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);


% create initial pressure distribution
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
%P0 = [M1 [M3; M2; M3] M4];
P0 = zeros(Nx,Ny); P0(floor(Nx*6/11), floor(7/8*Ny)) = 1;


source.p0 = P0;

% Define the sensors
nSensors = 3;
U = zeros(Nx, Ny);
U(floor(Nx/2), 25) = 1; % Sensor 1
U(floor(Nx*2/3), 25) = 1; % Sensor 2
U(25, floor(Ny/2)) = 1; % Sensor 3
sensor.mask = U; %ones(Nx, Ny); %U;

%%%%% Forward %%%%%%%%%%%%%%
% run the simulation
sensor.record = {'p', 'p_final'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% for visualisation, record pressure everywhere
sensor2 = sensor;
sensor2.mask = ones(Nx, Ny); %U;
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor2);

%%%%% Back propagation %%%%%
% reset the initial pressure
source.p0 = 0;
% assign the time reversal data
sensor.time_reversal_boundary_data = 1e1*sensor_data.p;
% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor);
% run the time-varying dirichlet source (I think this is the same as TR but allows me to record pressure time series
% everywhere)
source2 = source;
source2.p_mask = U;
source2.p_mode = 'dirichlet'; %'additive'; %'dirichlet'
source2.p = sensor.time_reversal_boundary_data(:,end:-1:1);
% %sensor.mask = ones(Nx, Ny);
% %sensor = rmfield(sensor,'time_reversal_boundary_data');
p0_recon2 = kspaceFirstOrder2D(kgrid, medium, source2, sensor2);
 
% Plot the forward propagating front overlayed over scaled speed of sound and source / sensor locations
figure, for j=1:10:411, imagesc(reshape(sensor_data2.p(:,j),Nx,Ny)/max(sensor_data2.p(:,j)) + U + P0 + medium.sound_speed/max(medium.sound_speed(:))); title(num2str(j)); pause; end
% Plot the time varying Dirichlet source 
figure, for j=1600:10:2164, imagesc(reshape(p0_recon2.p(:,j),Nx,Ny)/max(p0_recon2.p(:,j)) + U + P0 + medium.sound_speed/max(medium.sound_speed(:))); title(num2str(j)); pause; end

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
saveas(gcf, 'Example20_U_kWave.fig');


% plot the simulated sensor data
figure;
hold on;
for n = 1:nSensors
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMapV(n));
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, 'Example20_aSignal_kWave.fig');
