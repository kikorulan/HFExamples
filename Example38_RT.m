%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex38_gaussBeams;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex38_gaussBeams;

close all;
%clear all;

% Measure computational time
tic;
start_time = clock;

run colourMap;
load sensor_data.mat;
%========================================
% Grid definition
%========================================
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]

Rgrid = gridRT(Nx, dx, Ny, dy);
% Time signal
dt = 3e-7;
tMax = 4e-4;
Rgrid.setTime(dt, tMax);

% Sound speed
c0 = 1500;
Rgrid.setCMatrix(medium.sound_speed);

% Build initial pressure
Rgrid.setUMatrix(source.p0);

% Impulse Response
Rgrid.impulseResponse2D('IV', 1);

%========================================
% Ray Shooting
%========================================

% Number of rays & sources
nRays = 1000;
nSources = 3;

tMax = sqrt((Nx*dx)^2+(Ny*dy)^2)/c0;
tStep = dt;

% Sources locations
x1(1, 1, 1) = 2*Nx*dx/3;
x1(1, 1, 2) = 0; % Source 1
x2(1, 1, 1) = Nx*dx/4;
x2(1, 1, 2) = (Ny-1)*dy; % Source 2
x3(1, 1, 1) = 3*Nx*dx/4; 
x3(1, 1, 2) = (Ny-1)*dy; % Source 3

% Create the new sources
Rgrid.setDeltaX(Rgrid.dx/100);
Rgrid.newSource(x1, 0, pi, nRays, tStep, tMax);
Rgrid.newSource(x2, 0, -pi, nRays, tStep, tMax);
Rgrid.newSource(x3, 0, -pi, nRays, tStep, tMax);

% Sources
for n = 1:nSources
    Rgrid.computeHamil(n);
end

% Save results
save gridRT.mat Rgrid nRays nSources x1 x2 x3 -v7.3;

%==================================================================================
% Plot results
%==================================================================================
position = [700 700 450 600];
set(0,'DefaultFigurePaperPositionMode','auto');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex38_gaussBeams;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex38_gaussBeams;

%load gridRT.mat;
load sensor_data.mat;

%==============================
% Sound speed
%==============================
Rgrid.plot_soundSpeed();
title('');
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Sound Speed');
saveas(gcf, 'Example38_C', 'png'); 
saveas(gcf, 'Example38_C.fig'); 

%==============================
% Initial pressure
%==============================
figure;
flipGrid = Rgrid.u';
x = 0:Rgrid.dx:(Rgrid.Nx-1)*Rgrid.dx;
y = 0:Rgrid.dy:(Rgrid.Ny-1)*Rgrid.dy;
axis([0 x(end) 0 y(end)]);
hold on;
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Initial Pressure');
saveas(gcf, 'Example38_U', 'png'); 
saveas(gcf, 'Example38_U.fig'); 

%==============================
% Sensor Positions
%==============================
figure;
axis([0 x(end) 0 y(end)]);
hold on;
plot(x1(1), x1(2), 'or', 'MarkerSize', 10, 'LineWidth', 5);
plot(x2(1), x2(2), 'og', 'MarkerSize', 10, 'LineWidth', 5);
plot(x3(1), x3(2), 'ob', 'MarkerSize', 10, 'LineWidth', 5);
% Rays
for j = 1:20:nRays
    % Sources
    for n = 1:nSources
        plot(Rgrid.source(n).x(j, :, 1), Rgrid.source(n).x(j, :, 2), 'Color', colourMapV(n));
    end
end
leg = legend('Sensor 1', 'Sensor 2', 'Sensor 3');
set(leg, 'Position', [0.1 0.1 0.1 0.2])
xlabel('x (m)');
ylabel('y (m)');
title('Sensor positions');
set(gcf, 'pos', position);
saveas(gcf, 'Example38_SensorPositions.fig');

%==============================
% Isotime curve
%==============================
nCurves = 20;
xIso = Rgrid.findIsoTime(2, 20, 2e-4, 2.5e-4);
% Figure
figure;
flipGrid = Rgrid.u';
x = 0:Rgrid.dx:(Rgrid.Nx-1)*Rgrid.dx;
y = 0:Rgrid.dy:(Rgrid.Ny-1)*Rgrid.dy;
axis([0 x(end) 0 y(end)]);
hold on;
surf(x, y, flipGrid, 'EdgeColor', 'none');
hold on;
for n = 1:nCurves
    plot3(xIso(n, :, 1), xIso(n, :, 2), repmat(2, [1 size(xIso, 2)]), '-m');
end
for j = 1:20:nRays
    plot3(Rgrid.source(2).x(j, :, 1), Rgrid.source(2).x(j, :, 2), repmat(2, [1 Rgrid.source(2).nPoints]), '-g');
end
view(2);
legend('Isocurves for times tMin = 2e-4s, tMax = 2.5e-4s');
title('Isotime curve');
set(gcf, 'pos', position);
saveas(gcf, 'Example38_isotime', 'png');
saveas(gcf, 'Example38_isotime.fig');

%==============================
% Amplitude
%==============================
% Figure
Rgrid.plot_amplitude(2, 100, 'real');
title('Amplitude GB');
saveas(gcf, 'Example38_amplitude_GB', 'png');

%==============================
% Non-Filtered Signal
%==============================
% Time Signal Amplitude
figure;
hold on;
axis([0 3e-4 -.2 1.2]);
normBeam = max(Rgrid.source(1).aBeam);
for n = 1:nSources
    plot(Rgrid.source(n).tBeam, Rgrid.source(n).aBeam/normBeam, 'Color', colourMapV(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
xlabel('t (s)');
ylabel('Amplitude');
grid on;
box on;
title('Signal GB');
saveas(gcf, 'Example38_aSignal_GB', 'png');
saveas(gcf, 'Example38_aSignal_GB.fig');

%==============================
% Filtered Signal
%==============================
% Normalisation factor
normRT = max(Rgrid.source(1).aForward);
normKWave = max(sensor_data.p(1, :));

% Figure
figure;
%%% SIGNAL
axis([0 4e-4 -1.2 1.2]);
hold on;
for n = 1:nSources
    % Plot and compare
    plot(Rgrid.tForward, Rgrid.source(n).aForward/normRT, 'Color', colourMapV(n));
    plot(kgrid.t_array, sensor_data.p(n, :)/normKWave, 'Color', colourMapV(n+nSources), 'LineWidth', 2);
end
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
grid on;
box on;
title('Signal GB');
saveas(gcf, 'Example38_aSignalConv_GB', 'png');
saveas(gcf, 'Example38_aSignalConv_GB.fig');

%==============================
% Error
%==============================
figure;
axis([0 4e-4 -.2 .2]);
hold on;
for n = 1:nSources
    error = Rgrid.source(n).aForward/normRT - sensor_data.p(n, :)/normKWave;
    plot(Rgrid.tForward, error, 'Color', colourMapV(n));
end
xlabel('t (s)');
ylabel('Error');
legend('Error 1', 'Error 2', 'Error 3');
grid on;
box on;
title('Error GB');
saveas(gcf, 'Example38_error_GB', 'png');
saveas(gcf, 'Example38_error_GB.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/;
