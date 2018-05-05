% Monopole Point Source In A Homogeneous Propagation Medium Example
%

clear all;
close all;
run colorMap;
% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 60;           % number of grid points in the x (row) direction
Ny = 60;           % number of grid points in the y (column) direction
dx1 = 60e-3/Nx;    	% grid point spacing in the x direction [m]
dy1 = 60e-3/Ny;            % grid point spacing in the y direction [m]
dx2 = 30e-3/Nx;    	% grid point spacing in the x direction [m]
dy2 = 30e-3/Ny;            % grid point spacing in the y direction [m]
kgridImpulse1 = makeGrid(Nx, dx1, Ny, dy1);
kgridImpulse2 = makeGrid(Nx, dx2, Ny, dy2);
% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]


%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% define a single source point
source.p0 = zeros(Nx, Ny);
source.p0(end - Nx/4, Ny/2) = 1;

% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/4, Ny/2) = 1;

% define the acoustic parameters to record
sensor.record = {'p', 'p_final', 'u_final'};


nImpulses = 3;
% create the time array
for j = 1:nImpulses
    dt = 5e-8*10^(j/2-0.5)
    tArray{j} = 0:dt:1.2e-4;
    tArray{j+nImpulses} = 0:dt:1.2e-4;
    kgridImpulse1.t_array = tArray{j};
    kgridImpulse2.t_array = tArray{j+nImpulses};
    % run the simulation
    sensor_dataImpulse{j} = kspaceFirstOrder2D(kgridImpulse1, medium, source, sensor);
    sensor_dataImpulse{j+nImpulses} = kspaceFirstOrder2D(kgridImpulse2, medium, source, sensor);
end
% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulated sensor data
figure;
hold on;
for j = 1:nImpulses
    plot(tArray{j}, sensor_dataImpulse{j}.p, 'Color', colourMap(j));
end
xlabel('Time');
ylabel('Signal Amplitude');
grid on;
title('Sensor Pressure Signal');
legend('dt = 5e-8, dx = 1e-3', 'dt = 1.58e-7, dx = 1e-3', 'dt = 5e-7, dx = 1e-3', ...
 'dt = 5e-8, dx = 5e-4', 'dt = 1.58e-7, dx = 5e-4', 'dt = 5e-7, dx = 5e-4');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex15_kWave_imp/Example15_Impulse.fig'); 
