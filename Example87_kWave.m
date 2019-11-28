% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex87_lens;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
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
v1 = -0.3;
Z = kgrid.x.^2 + kgrid.y.^2;
sigma = 12*1e-4;
c = c0*(1+v1*exp(-Z/sigma/sigma));
medium.sound_speed = c;
medium.density = 1;

    
% compute time
dt = 1e-8;
tMax = 2e-5;
kgrid.t_array = 0:dt:tMax;
%[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% smooth the initial pressure distribution and restore the magnitude
sigmaU = 5e-4;
Z1 = (kgrid.x-30e-4).^2 + (kgrid.y+64e-4).^2;
source_low.p0 = exp(-Z1/sigmaU/sigmaU);
Z2 = (kgrid.x).^2 + (kgrid.y-25e-4).^2;
source_mid.p0 = exp(-Z2/sigmaU/sigmaU);
Z3 = (kgrid.x).^2 + (kgrid.y-77e-4).^2;
source_high.p0 = exp(-Z3/sigmaU/sigmaU);

%=========================================================================
% SIMULATION
%=========================================================================
% Define the sensors
sensor.mask = zeros(Nx, Ny);
sensor.mask(49, 1) = 1; % Sensor 1
sensor.mask(50, 1) = 1; % Sensor 2
sensor.mask(51, 1) = 1; % Sensor 3
sensor.mask(49, 2) = 1; % Sensor 4
sensor.mask(50, 2) = 1; % Sensor 5
sensor.mask(51, 2) = 1; % Sensor 6
sensor.mask(49, 3) = 1; % Sensor 7
sensor.mask(50, 3) = 1; % Sensor 8
sensor.mask(51, 3) = 1; % Sensor 9

sensor.record = {'p', 'u_non_staggered'};
% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PlotPML', false, 'Smooth', false};
% run the simulation
%sensor_data_low = kspaceFirstOrder2D(kgrid, medium, source_low, sensor, input_args{:});
%sensor_data_mid = kspaceFirstOrder2D(kgrid, medium, source_mid, sensor, input_args{:});
sensor_data_high = kspaceFirstOrder2D(kgrid, medium, source_high, sensor, input_args{:});

%save sensor_data.mat kgrid sensor source_low source_mid source_high medium c0 dt sensor_data_low sensor_data_mid sensor_data_high input_args;

%=========================================================================
% VISUALISATION
%=========================================================================

%==============================
% PRESSURE and SOUND SPEED
%==============================
figure;
surf(kgrid.y_vec, kgrid.x_vec, source_low.p0 +  source_mid.p0 + source_high.p0, 'EdgeColor', 'none');
view(2);
figure;
imagesc(medium.sound_speed)

%==============================
% Pressure
%==============================
figure;
plot(kgrid.t_array, sensor_data_low.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, sensor_data_mid.p(5, :), 'Color', 'g');
plot(kgrid.t_array, sensor_data_high.p(5, :), 'Color', 'b');
legend('IP low', 'IP mid', 'IP high');
xlabel('t (s)');
ylabel('amplitude');


%==============================
% Velocity
%==============================
figure;
plot(kgrid.t_array, sensor_data_low.ux_non_staggered(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, sensor_data_mid.ux_non_staggered(5, :), 'Color', 'g');
plot(kgrid.t_array, sensor_data_high.ux_non_staggered(5, :), 'Color', 'b');
legend('IP low', 'IP mid', 'IP high');
xlabel('t (s)');
ylabel('amplitude');

figure;
plot(kgrid.t_array, sensor_data_low.uy_non_staggered(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, sensor_data_mid.uy_non_staggered(5, :), 'Color', 'g');
plot(kgrid.t_array, sensor_data_high.uy_non_staggered(5, :), 'Color', 'b');
legend('IP low', 'IP mid', 'IP high');
xlabel('t (s)');
ylabel('amplitude');

%================================================================================
% VECTOR DECOMPOSITION
%================================================================================
value_matrix = @(sensor_data, index) [sensor_data(1, index) sensor_data(4, index) sensor_data(7, index); ...
                                      sensor_data(2, index) sensor_data(5, index) sensor_data(8, index); ...
                                      sensor_data(3, index) sensor_data(6, index) sensor_data(9, index)];
[X, Y] = meshgrid(-1:1);
Xq = @(angle) [sin(angle), -sin(angle), cos(angle), -cos(angle)];
Yq = @(angle) [cos(angle), -cos(angle), -sin(angle), sin(angle)];

sensor_data = sensor_data_low.p;
s1 = sensor_data(5, :);
angle_sequence = 0:pi/10:2*pi;
max_angle = 0*angle_sequence;

ii = 0;
for angle = angle_sequence
    ii = ii+1
    % Interpolation
    p_1 = 0*kgrid.t_array;
    p_0 = 0*kgrid.t_array;
    for index = 5:length(kgrid.t_array)-4
        Vq = interp2(X, Y, value_matrix(sensor_data, index), Xq(angle), Yq(angle));
        %Vq = [sensor_data(4, index) sensor_data(6, index) sensor_data(2, index) sensor_data(8, index)];
        p_vv  = (Vq(1) - 2*sensor_data(5, index) + Vq(2))/dx/dx;
        p_vvT = 0*(Vq(3) - 2*sensor_data(5, index) + Vq(4))/dx/dx;
        p_2 = p_vv + p_vvT;
        p_1(index) = p_1(index-1) + c0*c0*p_2*dt;
    end
    for index = 5:length(kgrid.t_array)-4
        p_0(index) = p_0(index-1) + p_1(index)*dt;
    end
    %p_0 = p_0/(abs(sin(angle)) + abs(cos(angle)));
    max_angle(ii) = max(p_0)/max(s1);
end
figure;
plot(angle_sequence, max_angle);
figure;
plot(kgrid.t_array, s1, 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');


%================================================================================
% DERIVE PRESSURE FROM VELOCITY
%================================================================================

% dtu - low
ux_x  = (sensor_data_low.ux_non_staggered(6, :) - sensor_data_low.ux_non_staggered(4, :))/dx/2;
uy_y  = (sensor_data_low.uy_non_staggered(8, :) - sensor_data_low.uy_non_staggered(2, :))/dx/2;
rho_t = -ux_x - uy_y; 
rho = 0*kgrid.t_array;
p_0 = 0*kgrid.t_array;
for ii = 3:length(kgrid.t_array)-2
    rho(ii) = rho(ii-1) + dt*(rho_t(ii)+rho_t(ii-1))/2; 
    p_0(ii) = c0*c0*rho(ii);
end
figure;
plot(kgrid.t_array, sensor_data_low.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');

% dtu - mid
ux_x  = (sensor_data_mid.ux_non_staggered(6, :) - sensor_data_mid.ux_non_staggered(4, :))/dx/2;
uy_y  = (sensor_data_mid.uy_non_staggered(8, :) - sensor_data_mid.uy_non_staggered(2, :))/dx/2;
rho_t = -ux_x - uy_y; 
rho = 0*kgrid.t_array;
p_0 = 0*kgrid.t_array;
for ii = 3:length(kgrid.t_array)-2
    rho(ii) = rho(ii-1) + dt*(rho_t(ii)+rho_t(ii-1))/2; 
    p_0(ii) = c0*c0*rho(ii);
end
figure;
plot(kgrid.t_array, sensor_data_mid.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');

% dtu - high
ux_x  = (sensor_data_high.ux_non_staggered(6, :) - sensor_data_high.ux_non_staggered(4, :))/dx/2;
uy_y  = (sensor_data_high.uy_non_staggered(8, :) - sensor_data_high.uy_non_staggered(2, :))/dx/2;
rho_t = -ux_x - uy_y; 
rho = 0*kgrid.t_array;
p_0 = 0*kgrid.t_array;
for ii = 3:length(kgrid.t_array)-2
    rho(ii) = rho(ii-1) + dt*(rho_t(ii)+rho_t(ii-1))/2; 
    p_0(ii) = c0*c0*rho(ii);
end

figure;
plot(kgrid.t_array, sensor_data_high.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');

%================================================================================
% ALTERNATIVE
%================================================================================
% dtu - low
ux_t_x1  = (sensor_data_low.ux_non_staggered(4, 3:end) - sensor_data_low.ux_non_staggered(4, 1:end-2))/dt/2;
ux_t_x2  = (sensor_data_low.ux_non_staggered(6, 3:end) - sensor_data_low.ux_non_staggered(6, 1:end-2))/dt/2;
ux_t_x = (ux_t_x2 - ux_t_x1)/dx/2;
uy_t_y1  = (sensor_data_low.uy_non_staggered(2, 3:end) - sensor_data_low.uy_non_staggered(2, 1:end-2))/dt/2;
uy_t_y2  = (sensor_data_low.uy_non_staggered(8, 3:end) - sensor_data_low.uy_non_staggered(8, 1:end-2))/dt/2;
uy_t_y = (uy_t_y2 - uy_t_y1)/dx/2;

p_2 = - ux_t_x - uy_t_y;
p_1 = 0*kgrid.t_array;
p_0 = 0*kgrid.t_array;
for ii = 5:length(kgrid.t_array)-4
    p_1(ii) = p_1(ii-1) + c0*c0*p_2(ii)*dt;
end
for ii = 5:length(kgrid.t_array)-4
    p_0(ii) = p_0(ii-1) + p_1(ii)*dt;
end
figure;
plot(kgrid.t_array, sensor_data_low.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');

% dtu - mid
ux_t_x1  = (sensor_data_mid.ux_non_staggered(4, 3:end) - sensor_data_mid.ux_non_staggered(4, 1:end-2))/dt/2;
ux_t_x2  = (sensor_data_mid.ux_non_staggered(6, 3:end) - sensor_data_mid.ux_non_staggered(6, 1:end-2))/dt/2;
ux_t_x = (ux_t_x2 - ux_t_x1)/dx/2;
uy_t_y1  = (sensor_data_mid.uy_non_staggered(2, 3:end) - sensor_data_mid.uy_non_staggered(2, 1:end-2))/dt/2;
uy_t_y2  = (sensor_data_mid.uy_non_staggered(8, 3:end) - sensor_data_mid.uy_non_staggered(8, 1:end-2))/dt/2;
uy_t_y = (uy_t_y2 - uy_t_y1)/dx/2;

p_2 = - ux_t_x - uy_t_y;
p_1 = 0*kgrid.t_array;
p_0 = 0*kgrid.t_array;
for ii = 5:length(kgrid.t_array)-4
    p_1(ii) = p_1(ii-1) + c0*c0*p_2(ii)*dt;
end
for ii = 5:length(kgrid.t_array)-4
    p_0(ii) = p_0(ii-1) + p_1(ii)*dt;
end
figure;
plot(kgrid.t_array, sensor_data_mid.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');

% dtu - high
ux_t_x1  = (sensor_data_high.ux_non_staggered(4, 3:end) - sensor_data_high.ux_non_staggered(4, 1:end-2))/dt/2;
ux_t_x2  = (sensor_data_high.ux_non_staggered(6, 3:end) - sensor_data_high.ux_non_staggered(6, 1:end-2))/dt/2;
ux_t_x = (ux_t_x2 - ux_t_x1)/dx/2;
uy_t_y1  = (sensor_data_high.uy_non_staggered(2, 3:end) - sensor_data_high.uy_non_staggered(2, 1:end-2))/dt/2;
uy_t_y2  = (sensor_data_high.uy_non_staggered(8, 3:end) - sensor_data_high.uy_non_staggered(8, 1:end-2))/dt/2;
uy_t_y = (uy_t_y2 - uy_t_y1)/dx/2;

p_2 = - ux_t_x - uy_t_y;
p_1 = 0*kgrid.t_array;
p_0 = 0*kgrid.t_array;
for ii = 5:length(kgrid.t_array)-4
    p_1(ii+1) = p_1(ii) + c0*c0*p_2(ii)*dt;
end
for ii = 5:length(kgrid.t_array)-4
    p_0(ii+1) = p_0(ii) + p_1(ii)*dt;
end
figure;
plot(kgrid.t_array, sensor_data_high.p(5, :), 'Color', 'r');
hold on;
plot(kgrid.t_array, p_0, 'Color', 'b');
