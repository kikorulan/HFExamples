% Impulse response for a medium with 1500 m/s sound speed

%close all;
%clear all;
load gridRT_Ray.mat;

%==========================================================================
% SIMULATION
%==========================================================================

% create the computational grid
Nx = 60;       % number of grid points in the x (row) direction      
Ny = 60;       % number of grid points in the y (column) direction   
dx = grid.dx;  % grid point spacing in the x direction [m]           
dy = grid.dy;  % grid point spacing in the y direction [m]           
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = grid.cFilter;  % [m/s]

% create the time array
tau = grid.source(1).step;
cMax = max(grid.c(:));
xMax = 2*Nx*dx;
tMax = xMax/medium.sound_speed;
dt = tau/cMax/cMax;
kgrid.t_array = 0:dt:tMax + dt;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);


% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, Ny/2) = 1;
% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};

% ===============================
% Initial Value Problem
% ===============================
% define a single source point
source_IV.p0 = zeros(Nx, Ny);
source_IV.p0(end - Nx/4, Ny/2) = 1;

% ===============================
% Run the simulation
% ===============================
sensorData_IV = kspaceFirstOrder2D(kgrid, medium, source_IV, sensor);

%================================
% OBTAIN THE FILTERS
%================================
% Initial Value source
maxF_IV = max(sensorData_IV.p);
delay_IV = find(sensorData_IV.p == maxF_IV) - 1;
flipFilter_IV = fliplr(sensorData_IV.p(1:delay_IV));
indexIni_IV = find(flipFilter_IV < 1e-3*maxF_IV, 1);
flipFilter_IV = fliplr(sensorData_IV.p(delay_IV+1:end));
indexFin_IV = find(abs(flipFilter_IV) < 1e-3*maxF_IV, 1);

Filter_IV = sensorData_IV.p(delay_IV-indexIni_IV+1:end-indexFin_IV);
tArray_IV = kgrid.t_array(1:end-indexFin_IV-delay_IV+indexIni_IV);


% =========================================================================
% VISUALISATION
% =========================================================================

% plot the Filters
figure;
[t_sc, scale, prefix] = scaleSI(max(tArray_IV(:)));
hold on;
plot(tArray_IV*scale, Filter_IV, 'b-');
xlabel(['Time [' prefix 's]']);
ylabel('Signal Amplitude');
axis tight;
title('Filter');
legend('IV');

%clear all;

% =========================================================================
% DELAY
% =========================================================================

% Theoretical delay
delay_Th = Nx/2*dx/medium.sound_speed/dt;
% Simulation delay
delay_kW = delay_IV;

% Difference
delay_diff = delay_kW - delay_Th;


%%  % ==============================================
%%  % OBTAIN THE FILTERS - CAUSAL IMPLEMENTATION 
%%  % ==============================================
%%  % Initial Value source
%%  maxF_IV = max(sensorData_IV.p);
%%  indexMax_IV = find(sensorData_IV.p == maxF_IV);
%%  flipFilter_IV = fliplr(sensorData_IV.p(1:indexMax_IV));
%%  indexMin_IV = find(flipFilter_IV < 0.01*maxF_IV, 1);
%%  Filter_IV = sensorData_IV.p(indexMax_IV-indexMin_IV+1:end);
%%  tArray_IV = kgrid.t_array(1:length(Filter_IV));
