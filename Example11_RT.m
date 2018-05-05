%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex11_source;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex11_source;

close all;
%clear all;

load kgrid_data.mat;

% Measure computational time
tic;
start_time = clock;
 
run colourMap;
%========================================
% Grid definition
%========================================
Nx = 240;           % number of grid points in the x (row) direction
Ny = 360;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]

Rgrid = gridRT(Nx, dx, Ny, dy);
% Sound speed
c0 = 1500;
% Auxiliary matrices
c = c0*ones(Nx, Ny);
Rgrid.setCMatrix(c);

% Initial pressure
Rgrid.setUMatrix(source.p0);
%========================================
% Impulse Response
%========================================
% Set time
dt = 3e-7;
tMax = 4e-4;
Rgrid.setTime(dt, tMax);

Rgrid.impulseResponse2D('IV');

save gridRT_impulse.mat Rgrid;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse.mat;
% Number of rays & sources
nRays = 1000;
nSources = 3;

tauMax = sqrt((Rgrid.Nx*Rgrid.dx)^2+(Rgrid.Ny*Rgrid.dy)^2)*c0;
tauStep = min(Rgrid.dx, Rgrid.dy)*c0/3;

tMax = sqrt((Rgrid.Nx*Rgrid.dx)^2+(Rgrid.Ny*Rgrid.dy)^2)/c0;
%tStep = min(Rgrid.dx, Rgrid.dy)/c0/3;
tStep = dt;

% Sources locations
x1(1, 1, 1) = Rgrid.Nx*Rgrid.dx/2; 
x1(1, 1, 2) = 0; 
x2(1, 1, 1) = 5*Rgrid.Nx*Rgrid.dx/6; 
x2(1, 1, 2) = 0; 
x3(1, 1, 1) = 0; 
x3(1, 1, 2) = Rgrid.Ny*Rgrid.dy/3; 

% Create the new sources
Rgrid.setDeltaX(Rgrid.dx/100);
Rgrid.newSource(x1, 0, pi, nRays, tStep, tMax);
Rgrid.newSource(x2, 0, pi, nRays, tStep, tMax);
Rgrid.newSource(x3, -pi/2, pi/2, nRays, tStep, tMax);
  
% Sources
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    Rgrid.computeHamil(n);
end

% Save results
save gridRT_ray.mat Rgrid nRays nSources x1 x2 x3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex11_source;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex11_source;

%load gridRT.mat;
load kgrid_data.mat;
  
x = 0:Rgrid.dx:(Rgrid.Nx-1)*Rgrid.dx;
y = 0:Rgrid.dy:(Rgrid.Ny-1)*Rgrid.dy;

position = [700 700 480 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%%  %==============================
%%  % Sound speed
%%  %==============================
%%  Rgrid.plotSoundSpeed();
%%  saveas(gcf, 'Example11_C', 'epsc'); 
%%  saveas(gcf, 'Example11_C.fig'); 

%==============================
% Initial pressure
%==============================
figure;
flipGrid = Rgrid.u';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
axis([0 Rgrid.Nx*Rgrid.dx 0 Rgrid.Ny*Rgrid.dy]);
xlabel('x (m)');
ylabel('y (m)');
colorbar();

set(gcf, 'pos', position);
saveas(gcf, 'Example11_U', 'png');
saveas(gcf, 'Example11_U.fig'); 

%==============================
% Sensor positions
%==============================
figure;
axis([0 x(end) 0 y(end)]);
hold on;
plot(x1(1), x1(2), 'or', 'MarkerSize', 10, 'LineWidth', 6);
plot(x2(1), x2(2), 'og', 'MarkerSize', 10, 'LineWidth', 6);
plot(x3(1), x3(2), 'ob', 'MarkerSize', 10, 'LineWidth', 6);
xlabel('x (m)');
ylabel('y (m)');
grid on;
box on;
set(gcf, 'pos', [100 100 400 600]);
leg = legend('Sensor 1', 'Sensor 2', 'Sensor 3');
set(leg,'Position',[0.5 0.5 0.1 0.2]);
saveas(gcf, 'Example11_SensorPositions', 'png');
saveas(gcf, 'Example11_SensorPositions.fig'); 


%%  %==============================
%%  % Sources and their respective Rays
%%  %==============================
%%  figure;
%%  hold on;
%%  axis([0 Rgrid.Nx*Rgrid.dx 0 Rgrid.Ny*Rgrid.dy]);
%%  % Sources
%%  for n = 1:nSources
%%      % Rays
%%      for j = 1:10:nRays
%%          plot(Rgrid.source(n).x(j, :, 1), Rgrid.source(n).x(j, :, 2), 'Color', colourMapV(n));
%%      end
%%  end
%%  legend('Sensor 1', 'Sensor 2', 'Sensor 3');
%%  saveas(gcf, 'Example11_Ray', 'epsc');
%%  saveas(gcf, 'Example11_Ray.fig');

%==============================
% Non-Filtered Signal
%==============================
normBeam = max(Rgrid.source(1).aBeam);
% Time Signal Amplitude
figure;
hold on;
for n = 1:nSources
    plot(Rgrid.source(n).tBeam, Rgrid.source(n).aBeam/normBeam, 'Color', colourMapV(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
xlabel('t (s)');
ylabel('Amplitude');
grid on;
box on;
saveas(gcf, 'Example11_aSignal', 'epsc');
saveas(gcf, 'Example11_aSignal.fig');

%==============================
% Filtered Signal
%==============================
% Normalisation factor
normRT = max(Rgrid.source(1).aForward);
normKWave = max(sensor_data.p(1, :));

% Figure
figure;
hold on;
axis([0 4e-4 -1.2 1.2]);
for n = 1:nSources
    % Plot and compare
    plot(Rgrid.tForward, Rgrid.source(n).aForward/normRT, 'Color', colourMapV(n));
    plot(kgrid.t_array, sensor_data.p(n, :)/normKWave, 'Color', colourMapV(n+nSources), 'LineWidth', 2);
end
grid on;
box on;
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
saveas(gcf, 'Example11_aSignalConv', 'epsc');
saveas(gcf, 'Example11_aSignalConv.fig');
 
%%  %==============================
%%  % Non-Filtered + Filtered Signal
%%  %=============================
%%  normBeam = max(Rgrid.source(1).aBeam);
%%  normRT = max(Rgrid.source(1).aForward);
%%  % Figure
%%  figure;
%%  hold on;
%%  for n = 1:nSources
%%      % Plot and compare
%%      plot(Rgrid.tForward, Rgrid.source(n).aForward/normRT, 'Color', colourMapV(n), 'LineWidth', 2);
%%      plot(Rgrid.source(n).tBeam, Rgrid.source(n).aBeam/normBeam, 'Color', colourMapV(n+nSources));
%%  end
%%  
%%  legend('Sensor 1 RT', 'Sensor 1 kWave', ...
%%      'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
%%  saveas(gcf, 'Example11_aSignal_both.fig');

%==============================
% Error
%=============================
figure;
hold on;
axis([0 4e-4 -.2 .2]);
for n = 1:nSources
    % Plot and compare
    errorRT{n} = Rgrid.source(n).aForward/normRT - sensor_data.p(n, :)/normKWave;
    plot(Rgrid.tForward, errorRT{n}, 'Color', colourMapV(n));
end
grid on;
box on;
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, 'Example11_error', 'epsc');
saveas(gcf, 'Example11_error.fig');


%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/;


