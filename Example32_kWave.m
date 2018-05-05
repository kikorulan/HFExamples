% Homogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex32_ampl;

clear all;
close all;

run colourMap;
%=========================================================================
% SIMULATION
%=========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
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
dt = 3e-8;
tMax = 1.6e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% Build initial pressure
u0 = zeros(Nx, Ny);
u0(64, 32) = 1;
source.p0 = u0;
% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(kgrid, source.p0, true);

%sensor.mask = gridAux.u;
sensor.mask = zeros(Nx, Ny);
sensor.mask(96, 40:8:224) = 1;
% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
nSensors = size(sensor_data, 1);

save sensor_data.mat sensor_data kgrid source nSensors;

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex32_ampl;

% plot the simulated sensor data
nColours = cool(nSensors);
figure;
hold on;
for n = 1:nSensors;
    plot(kgrid.t_array, sensor_data(n, :), 'Color', nColours(n, :));
end
ylabel('Sensor');
xlabel('Time Step');
saveas(gcf, 'Example32_aSignal_kWave.fig');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex32_ampl;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
