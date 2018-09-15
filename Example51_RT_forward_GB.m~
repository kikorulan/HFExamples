%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex51_reconstruction2D;

close all;
%clear all;

load sensor_data.mat;

run colourMap;
%========================================
% Rgrid definition
%========================================
Nx = 128;           % number of Rgrid points in the x (row) direction
Ny = 256;           % number of Rgrid points in the y (column) direction
dx = 2e-4;        % Rgrid point spacing in the x direction [m]
dy = 2e-4;        % Rgrid point spacing in the y direction [m]
Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
Rgrid.setCMatrix(medium.sound_speed);

% Set initial pressure
Rgrid.setUMatrix(sourceKW.p0);

%========================================
% Impulse Response
%========================================
% Set time
dt = 5e-8;
tMax = 4e-5;
Rgrid.setTime(dt, tMax);
% Compute impulse response
Rgrid.impulse_additive('IV');

save gridRT_impulse.mat Rgrid;

%========================================
% Ray Shooting Parameters
%========================================
load gridRT_impulse;
cMax = max(Rgrid.c(:));
dt = 5e-8;
%dt = min(Rgrid.dx, Rgrid.dy)/cMax/2;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 2000;% 800
nSources = 764;%256

% Parametrisation
tMax = sqrt((Rgrid.Nx*Rgrid.dx)^2+(Rgrid.Ny*Rgrid.dy)^2)/cMax;
tStep = dt;

%========================================
% Forward Problem
%========================================
% Sources locations
clear x;
for n = 1:Rgrid.Nx
    x{n} = cat(3, (n-1)*Rgrid.dx, 0);
end
for n = 1:Rgrid.Ny-2
    x{2*n-1+Rgrid.Nx}   = cat(3,                   0, n*Rgrid.dy);
    x{2*n  +Rgrid.Nx}   = cat(3, (Rgrid.Nx-1)*Rgrid.dx, n*Rgrid.dy);
end
for n = 1:Rgrid.Nx
    x{n+Rgrid.Nx+2*(Rgrid.Ny-2)} = cat(3, (n-1)*Rgrid.dx, (Rgrid.Ny-1)*Rgrid.dy);
end
source = Rgrid.computeForwardParallel(x, 0, 2*pi-0.01, nRays, tStep, tMax, true);
aForward_RT = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    aForward_RT(n, :) = source(n).aForward;
end

%%  %========================================
%%  % Sensor Selection
%%  %========================================
%%  % Sources locations
%%  n1 = round(Rgrid.Nx + 2*(round(Rgrid.Ny/3) - 2)   + 2);     % 1st sensor: Nx + 2*(Ny/3-2)   + 2
%%  n2 = round(Rgrid.Nx + 2*(round(2*Rgrid.Ny/3) - 2) + 1);     % 2nd sensor: Nx + 2*(2*Ny/3-2) + 1
%%  n3 = round(Rgrid.Nx + 2*(Rgrid.Ny - 2) + round(Rgrid.Nx/2)); % 3rd sensor: Nx + 2*(Ny-2)     + Nx/2
%%  sensor1 = x{n1};
%%  sensor2 = x{n2};
%%  sensor3 = x{n3};
%%  sourceSel(1) = Rgrid.newSource(sensor1, pi/2, 3*pi/2, nRays, tStep, tMax);
%%  sourceSel(2) = Rgrid.newSource(sensor2, -pi/2, pi/2, nRays, tStep, tMax);
%%  sourceSel(3) = Rgrid.newSource(sensor3, pi, 2*pi, nRays, tStep, tMax);
%%  Rgrid.computeHamil(sourceSel(1));
%%  Rgrid.computeHamil(sourceSel(2));
%%  Rgrid.computeHamil(sourceSel(3));

% Save results
%save gridRT.mat Rgrid nRays source x -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex51_reconstruction2D;

%load gridRT.mat;
%load sensor_data.mat;

%%  axisRgrid = [0 1e3*Rgrid.xAxis(end) 0 1e3*Rgrid.yAxis(end)];
%%  position       = [700 700 320 630];
%%  positionNoY    = [700 700 300 600];
%%  positionNoYBar = [700 700 363 600];
%%  positionYBar   = [700 700 390 630];
%%  
%%  %==============================
%%  % Initial Pressure
%%  %==============================
%%  figure;
%%  surf(1e3*Rgrid.xAxis, 1e3*Rgrid.yAxis, Rgrid.u', 'EdgeColor', 'none');
%%  view(2);
%%  saveas(gcf, 'Example51_U.fig'); 
%%  
%%  %==============================
%%  % Sources and their respective Rays
%%  %==============================
%%  nRaysPlot = 140;
%%  colorList = cool(nRaysPlot);
%%  % Sensor 1
%%  figure;
%%  hold on;
%%  axis(axisGrid);
%%  for n = 1:nRaysPlot
%%      index = n*floor(nRays/nRaysPlot);
%%      plot(1e3*sourceSel(1).x(index, :, 1), 1e3*sourceSel(1).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
%%  end
%%  plot(1e3*sensor1(1), 1e3*sensor1(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','r');
%%  box on;
%%  set(gcf, 'pos', position);
%%  xlabel('x [mm]');
%%  ylabel('y [mm]');
%%  saveas(gcf, 'Example51_Sensor1_rays', 'png');
%%  saveas(gcf, 'Example51_Sensor1_rays.fig'); 
%%  
%%  % Sensor 2
%%  figure;
%%  hold on;
%%  axis(axisGrid);
%%  for n = 1:nRaysPlot
%%      index = n*floor(nRays/nRaysPlot);
%%      plot(1e3*sourceSel(2).x(index, :, 1), 1e3*sourceSel(2).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
%%  end
%%  plot(1e3*sensor2(1), 1e3*sensor2(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','g');
%%  box on;
%%  set(gcf, 'pos', position);
%%  xlabel('x [mm]');
%%  saveas(gcf, 'Example51_Sensor2_rays', 'png');
%%  saveas(gcf, 'Example51_Sensor2_rays.fig'); 
%%  
%%  
%%  % Sensor 3
%%  figure;
%%  hold on;
%%  axis(axisGrid);
%%  for n = 1:nRaysPlot
%%      index = n*floor(nRays/nRaysPlot);
%%      plot(1e3*sourceSel(3).x(index, :, 1), 1e3*sourceSel(3).x(index, :, 2), 'Color', colorList(n, :), 'LineWidth', 1.5);
%%  end
%%  plot(1e3*sensor3(1), 1e3*sensor3(2), 'ok', 'MarkerSize', 8, 'LineWidth', 1, 'MarkerFaceColor','b');
%%  box on;
%%  set(gcf, 'pos', position);
%%  xlabel('x [mm]');
%%  saveas(gcf, 'Example51_Sensor3_rays', 'png');
%%  saveas(gcf, 'Example51_Sensor3_rays.fig'); 

%%  %==============================
%%  % Beam Signals
%%  %==============================
%%  figure;
%%  hold on;
%%  for n = 1:nSources
%%      plot(Rgrid.source(n).tBeam, Rgrid.source(n).aBeam);
%%  end
  
%%  %==============================
%%  % Time Signals
%%  %==============================
%%  normRT = max(sourceSel(2).aForward);
%%  normKWave = max(sensor_data(n2, :));
%%  % Signals
%%  signalRT1 = sourceSel(1).aForward/normRT;
%%  signalRT2 = sourceSel(2).aForward/normRT;
%%  signalRT3 = sourceSel(3).aForward/normRT;
%%  signalKW1 = sensor_data(n1, :)/normKWave;
%%  signalKW2 = sensor_data(n2, :)/normKWave;
%%  signalKW3 = sensor_data(n3, :)/normKWave;
%%  
%%  % Sensor 1
%%  figure; hold on;
%%  axis([0 40 -.6 1.2]);
%%  plot(1e6*Rgrid.tForward, signalRT1, 'Color', colourMapV(1), 'LineWidth', 3);
%%  plot(1e6*kgrid.t_array,  signalKW1, 'Color',           'k', 'LineWidth', 1.5);
%%  box on; grid on;
%%  legend('Sensor 1 - RT', 'Sensor 1 - kWave');
%%  xlabel('t [\mus]');
%%  ylabel('Amplitude');
%%  saveas(gcf, 'Example51_Sensor1_signal', 'png');
%%  saveas(gcf, 'Example51_Sensor1_signal.fig'); 
%%  % Sensor 2;
%%  figure; hold on;
%%  axis([0 40 -.6 1.2]);
%%  plot(1e6*Rgrid.tForward, signalRT2, 'Color', colourMapV(2), 'LineWidth', 3);
%%  plot(1e6*kgrid.t_array,  signalKW2, 'Color',           'k', 'LineWidth', 1.5);
%%  box on; grid on;
%%  legend('Sensor 2 - RT', 'Sensor 2 - kWave');
%%  xlabel('t [\mus]');
%%  ylabel('Amplitude');
%%  saveas(gcf, 'Example51_Sensor2_signal', 'png');
%%  saveas(gcf, 'Example51_Sensor2_signal.fig'); 
%%  % Sensor 3
%%  figure; hold on;
%%  axis([0 40 -.6 1.2]);
%%  plot(1e6*Rgrid.tForward, signalRT3, 'Color', colourMapV(3), 'LineWidth', 3);
%%  plot(1e6*kgrid.t_array,  signalKW3, 'Color',           'k', 'LineWidth', 1.5);
%%  box on; grid on;
%%  legend('Sensor 3 - RT', 'Sensor 3 - kWave');
%%  xlabel('t [\mus]');
%%  ylabel('Amplitude');
%%  saveas(gcf, 'Example51_Sensor3_signal', 'png');
%%  saveas(gcf, 'Example51_Sensor3_signal.fig'); 
%%  % Errors
%%  error1 = signalRT1 - signalKW1;
%%  error2 = signalRT2 - signalKW2;
%%  error3 = signalRT3 - signalKW3;
%%  figure; hold on;
%%  axis([0 40 -.6 1.2]);
%%  plot(1e6*Rgrid.tForward, error1, 'Color', colourMapV(1), 'LineWidth', 2);
%%  plot(1e6*Rgrid.tForward, error2, 'Color', colourMapV(2), 'LineWidth', 2);
%%  plot(1e6*Rgrid.tForward, error3, 'Color', colourMapV(3), 'LineWidth', 2);
%%  box on; grid on;
%%  legend('Error sensor 1', 'Error sensor 2', 'Error sensor 3');
%%  xlabel('t [\mus]');
%%  ylabel('Amplitude');
%%  saveas(gcf, 'Example51_error_signal', 'png');
%%  saveas(gcf, 'Example51_error_signal.fig'); 

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save RgridRT.mat Rgrid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

