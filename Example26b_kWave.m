% Find how is impulse response affected by change of medium

% =========================================================================
% SIMULATION
% =========================================================================
close all;
% create the computational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 0.1/Nx;        % grid point spacing in the x direction [m]
dy = 0.15/Ny;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

%%%%%%%%%%%%%%%%%%
% define the properties of the propagation medium
%%%%%%%%%%%%%%%%%
% Build domain
c0 = 1500;
v1 = -0.5;
v2 = 0.5;
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


%===================================
% Sensors
%===================================
% Sensor mask
sensor.mask = zeros(Nx, Ny);
x1 = 60;
x2 = 180;
y1 = 60;
y2 = 270;

sensor.mask(x1, y2) = 1;
sensor.mask(x2, y1) = 1;
sensor.mask(x2, y2) = 1;
% run the simulation
sensor.record = {'p', 'p_final'};

%===================================
% Build the source
%===================================
source.p0 = zeros(Nx, Ny);
source.p0(121, 73) = 1;

% compute time
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
%% Run the simulation
tArrayDef = kgrid.t_array;
%sensorDataDef = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% compute time
dt = 1e-8;
kgrid.t_array = 0:dt:3e-5;
tArray1 = kgrid.t_array;
%% Run the simulation
%sensorData1 = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% compute time
dt = 5e-9;
kgrid.t_array = 0:dt:3e-5;
tArray2 = kgrid.t_array;
%% Run the simulation
%sensorData2 = kspaceFirstOrder2D(kgrid, medium, source, sensor);


%==========================================================================
% PLOT RESULTS
%==========================================================================
x = 0:dx:(Nx-1)*dx;
y = 0:dy:(Ny-1)*dy;

% Plot sound speed
figure;
flipGrid = medium.sound_speed';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);

% Sensors
figure;
flipGrid = sensor.mask';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);

% Signals at Sensors
figure;
%axisY = [0 3e-5 -0.4 .8];
% RT Signals
%axis(axisY);
hold on;
plot(tArrayDef, sensorDataDef.p(1, :), 'Color', [0 0.6 0]);  
plot(tArrayDef, sensorDataDef.p(2, :), 'b');                 
plot(tArrayDef, sensorDataDef.p(3, :), 'Color', [0.7 0.7 0]);
legend('x21 - kW', 'x12 - kW', 'x22 - kW');
title('kWave signals - Default dt = 5.55e-8');
grid on;
%%  % RT Signals
%%  subplot(2, 1, 2);
%%  axis(axisY);
%%  hold on;
%%  plot(tArray1, sensorData1.p(1, :), 'Color', [0 0.6 0]);  
%%  plot(tArray1, sensorData1.p(2, :), 'b');                 
%%  plot(tArray1, sensorData1.p(3, :), 'Color', [0.7 0.7 0]);
%%  xlabel('t (s)');
%%  legend('x21 - kW', 'x12 - kW', 'x22 - kW');
%%  title('kWave signals - dt = 1e-8');
saveas(gcf, 'Example26b_kWave_dt.fig');
saveas(gcf, 'Example26b_kWave_dt', 'epsc');

%%  axis(axisY);
%%  hold on;
%%  plot(tArray1, sensorData2.p(1, :), 'Color', [0 0.6 0]);  
%%  plot(tArray1, sensorData2.p(2, :), 'b');                 
%%  plot(tArray1, sensorData2.p(3, :), 'Color', [0.7 0.7 0]);
%%  xlabel('t (s)');
%%  legend('x21 - kW', 'x12 - kW', 'x22 - kW');
%%  title('kWave signals - dt = 1e-9');
%%  saveas(gcf, 'Example26b_kWave_dt.fig');
%%  saveas(gcf, 'Example26b_kWave_dt', 'epsc');
