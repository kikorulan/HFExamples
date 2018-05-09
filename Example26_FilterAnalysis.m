close all;

%load Filter;
%load gridKWave.mat;
%load grid_IV.mat;

x1 = 60;
x2 = 180;
y1 = 60;
y2 = 270;

%==============================
% RT Signals
%==============================
signalRT11 = permute(grid_IV.pixelAmplitudeSignal(x1, y1, :), [1 3 2]);
signalRT21 = permute(grid_IV.pixelAmplitudeSignal(x2, y1, :), [1 3 2]);
signalRT12 = permute(grid_IV.pixelAmplitudeSignal(x1, y2, :), [1 3 2]);
signalRT22 = permute(grid_IV.pixelAmplitudeSignal(x2, y2, :), [1 3 2]);

% Norms
L2RT11 = sqrt(sum(signalRT11.*signalRT11));
L2RT12 = sqrt(sum(signalRT12.*signalRT12));
L2RT21 = sqrt(sum(signalRT21.*signalRT21));
L2RT22 = sqrt(sum(signalRT22.*signalRT22));
%==============================
% kWave Signals
%==============================
signalKW11 = permute(gridKWave.signal(x1, y1, :), [1 3 2]);
signalKW21 = permute(gridKWave.signal(x2, y1, :), [1 3 2]);
signalKW12 = permute(gridKWave.signal(x1, y2, :), [1 3 2]);
signalKW22 = permute(gridKWave.signal(x2, y2, :), [1 3 2]);

% Norms
L2KW11 = sqrt(sum(signalKW11.*signalKW11));
L2KW12 = sqrt(sum(signalKW12.*signalKW12));
L2KW21 = sqrt(sum(signalKW21.*signalKW21));
L2KW22 = sqrt(sum(signalKW22.*signalKW22));
%==============================
% Plots
%==============================
timeSignalRT = grid_IV.timeSignal;
timeSignalKW = gridKWave.time;

factorRT = L2KW21;
factorKW = L2RT21;

figure;
axisY = [0 1e-4 -0.4 .8];
% RT Signals
subplot(2, 1, 1);
axis(axisY);
set(gca, 'ytick', -0.5:0.1:1)
hold on;
%plot(timeSignalRT, factorRT*signalRT11, 'r');
plot(timeSignalRT, factorRT*signalRT21, 'Color', [0 0.6 0]); 
plot(timeSignalRT, factorRT*signalRT12, 'b');
plot(timeSignalRT, factorRT*signalRT22, 'Color', [0.7 0.7 0]); % dark yellow
grid on;
legend('x21 - RT', 'x12 - RT', 'x22 - RT');
xlabel('t (s)');
title('RT signals');
% kWave Signals
subplot(2, 1, 2);
axis(axisY);
set(gca, 'ytick', -0.5:0.1:1)
hold on;
%plot(timeSignal, factorKW*signalKW11, 'LineWidth', 1, 'Color', 'r');
plot(timeSignalKW, factorKW*signalKW21, 'LineWidth', 1, 'Color', [0 0.6 0]); % green
plot(timeSignalKW, factorKW*signalKW12, 'LineWidth', 1, 'Color', [0 0 1]);
plot(timeSignalKW, factorKW*signalKW22, 'LineWidth', 1, 'Color', [0.7 0.7 0]);
grid on;
legend('x21 - kW', 'x12 - kW', 'x22 - kW');
xlabel('t (s)');
title('kWave signals');
saveas(gcf, 'Example26_Signals.fig');
saveas(gcf, 'Example26_Signals', 'epsc'); 

%==============================
% Sensor Matrix
%==============================
x = 0:grid_IV.dx:(grid_IV.Nx-1)*grid_IV.dx;
y = 0:grid_IV.dy:(grid_IV.Ny-1)*grid_IV.dy;
sensorMatrix = zeros(grid_IV.Nx, grid_IV.Ny);
sensorMatrix(3:end-2, 3:end-2) = inf;
sensorMatrix(x2-2:x2+2, y1-2:y1+2) = 2;
sensorMatrix(x1-2:x1+2, y2-2:y2+2) = 1;
sensorMatrix(x2-2:x2+2, y2-2:y2+2) = 3;
sensorMatrix(119:123, 72:75) = 0;
% Border
figure;
flipGrid = sensorMatrix';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
title('Sensors');
saveas(gcf, 'Example26_Sensors.fig');
saveas(gcf, 'Example26_Sensors', 'epsc');

%==============================
% Sound speed
%==============================
grid_IV.plotSoundSpeed;
title('Sound Speed');
xlabel('x (m)');
ylabel('y (m)');
saveas(gcf, 'Example26_C.fig');

%==============================
% Attenuation
%==============================
attRT11 = grid_IV.pixelAttenuation(x1, y1);
msg = strcat({'Attenuation x11 = '}, num2str(attRT11));
disp(msg{1});
attRT21 = grid_IV.pixelAttenuation(x2, y1);
msg = strcat({'Attenuation x21 = '}, num2str(attRT21));
disp(msg{1});
attRT12 = grid_IV.pixelAttenuation(x1, y2);
msg = strcat({'Attenuation x12 = '}, num2str(attRT12));
disp(msg{1});
attRT22 = grid_IV.pixelAttenuation(x2, y2);
msg = strcat({'Attenuation x22 = '}, num2str(attRT22));
disp(msg{1});

%==============================
% Filter analysis
%==============================
figure;
plot(grid_IV.tFilter, grid_IV.Filter);
title('Filter')
saveas(gcf, 'Example26_Filter.fig');

%%  % Generate Filter time array
%%  % Set filter characteristics
%%  grid_IV.setCFilter(c0);
%%  grid_IV.setTFilter(tArray);
%%  grid_IV.setFilter(Filter);
%%  % Minimum and maximum sound speed
%%  cMin = 1200;
%%  cMax = 1800;
%%  tMin = grid_IV.source(1).tSignal(2); % time signal increment
%%  aSignal = grid_IV.source(1).aSignal;
%%  lengthASignal = length(aSignal);
%%  tFilter = 0:tMin:grid_IV.tFilter(end)*grid_IV.cFilter/cMax+tMin; % time array for the filter
%%  lengthFilter = length(tFilter);
%%  nFilters = 100; % number of filters;
%%  c = cMin:(cMax - cMin)/nFilters:cMax;
%%  FilterFFTMatrix = zeros(grid_IV.Nx, grid_IV.Ny, lengthFilter + lengthASignal);
%%  cMatrix = repmat(grid_IV.c, [1 1 lengthFilter + lengthASignal]);
%%  for i = 1:length(c)
%%      disp(int2str(i));
%%      % Compute splines for filter
%%      splineFilter = spline(grid_IV.tFilter*grid_IV.cFilter/c(i), grid_IV.Filter, tFilter);
%%      FilterPad = padarray(splineFilter, [0 lengthASignal], 0, 'post');
%%      FilterFFTaux = repmat(permute(fft(FilterPad), [1 3 2]), [grid_IV.Nx gridy.Ny 1]);
%%      binaryC = (cMatrix >= c(i));
%%      FilterFFTMatrix = FilterFFTMatrix.*(~binaryC) + FilterFFTaux.*binaryC;
%%  end
%%  FilterMatrix = ifft(FilterFFTMatrix, lengthFilter + lengthASignal, 3);
%%  
%%  Filter11 = permute(FilterMatrix(x1, y1, :), [1 3 2]);
%%  Filter21 = permute(FilterMatrix(x2, y1, :), [1 3 2]);
%%  Filter12 = permute(FilterMatrix(x1, y2, :), [1 3 2]);
%%  Filter22 = permute(FilterMatrix(x2, y2, :), [1 3 2]);
%%  
%%  figure;
%%  hold on;
%%  plot(Filter11);
%%  plot(Filter21);
%%  plot(Filter12);
%%  plot(Filter22);
%%  legend('x11', 'x21', 'x12', 'x22');
%%  
