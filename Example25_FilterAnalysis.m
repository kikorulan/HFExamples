close all;

load Filter;
%load gridKWave_TS.mat;
%%  load gridKWave_IV.mat;
%load grid_TS.mat;
%%  load grid_IV.mat;

x1 = 60;
x2 = 180;
y1 = 60;
y2 = 180;

%==============================
% RT Signals
%==============================
signalRT11 = permute(grid_TS.pixelAmplitudeSignal(x1, y1, :), [1 3 2]);
signalRT21 = permute(grid_TS.pixelAmplitudeSignal(x2, y1, :), [1 3 2]);
signalRT12 = permute(grid_TS.pixelAmplitudeSignal(x1, y2, :), [1 3 2]);
signalRT22 = permute(grid_TS.pixelAmplitudeSignal(x2, y2, :), [1 3 2]);

% Norms
L2RT11 = sqrt(sum(signalRT11.*signalRT11));
L2RT12 = sqrt(sum(signalRT12.*signalRT12));
L2RT21 = sqrt(sum(signalRT21.*signalRT21));
L2RT22 = sqrt(sum(signalRT22.*signalRT22));
%==============================
% kWave Signals
%==============================
signalKW11 = permute(gridKWave_TS(x1, y1, :), [1 3 2]);
signalKW21 = permute(gridKWave_TS(x2, y1, :), [1 3 2]);
signalKW12 = permute(gridKWave_TS(x1, y2, :), [1 3 2]);
signalKW22 = permute(gridKWave_TS(x2, y2, :), [1 3 2]);

% Norms
L2KW11 = sqrt(sum(signalKW11.*signalKW11));
L2KW12 = sqrt(sum(signalKW12.*signalKW12));
L2KW21 = sqrt(sum(signalKW21.*signalKW21));
L2KW22 = sqrt(sum(signalKW22.*signalKW22));
%==============================
% Plots
%==============================
timeSignal = grid_TS.timeSignal;

figure;
hold on;
factorRT = L2KW11;
factorKW = L2RT11;
% RT Signals
%plot(timeSignal, factorRT*signalRT11, 'r');
plot(timeSignal, factorRT*signalRT21, 'b');
plot(timeSignal, factorRT*signalRT12, 'g');
plot(timeSignal, factorRT*signalRT22, 'm');
% kWave Signals
%plot(timeSignal, factorKW*signalKW11, 'LineWidth', 2.5, 'Color', 'r');
plot(timeSignal, factorKW*signalKW21, 'LineWidth', 2.5, 'Color', 'b');
plot(timeSignal, factorKW*signalKW12, 'LineWidth', 2.5, 'Color', 'g');
plot(timeSignal, factorKW*signalKW22, 'LineWidth', 2.5, 'Color', 'm');
legend('x21 - RT', 'x12 - RT', 'x22 - RT', ...
         'x21 - kW', 'x12 - kW', 'x22 - kW');
title('RT vs kWave signals');
saveas(gcf, 'Example25_Signals.fig'); 

%==============================
% Sensor Matrix
%==============================
x = 1:1:grid_TS.Nx;
y = 1:1:grid_TS.Ny;
sensorMatrix = zeros(grid_TS.Nx, grid_TS.Ny);
sensorMatrix(x2, y1) = 1;
sensorMatrix(x1, y2) = 1;
sensorMatrix(x2, y2) = 1;
sensorMatrix(121, 49) = 2;
figure;
flipGrid = sensorMatrix';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
title('Sensors');
saveas(gcf, 'Example25_Sensors.fig'); 

%==============================
% Signal
%==============================
figure;
sineSignal = sin(2*pi./tArray(end).*tArray);
plot(tArray, sineSignal);
title('Source');
saveas(gcf, 'Example25_Source.fig');

%==============================
% Sound speed
%==============================
x = 1:1:grid_TS.Nx;
y = 1:1:grid_TS.Ny;
figure;
flipGrid = grid_TS.c';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
colorbar();
title('Sound Speed');
saveas(gcf, 'Example25_C.fig'); 

%==============================
% Attenuation
%==============================
attRT11 = grid_TS.pixelAttenuation(x1, y1);
msg = strcat({'Attenuation x11 = '}, num2str(attRT11));
disp(msg{1});
attRT21 = grid_TS.pixelAttenuation(x2, y1);
msg = strcat({'Attenuation x21 = '}, num2str(attRT21));
disp(msg{1});
attRT12 = grid_TS.pixelAttenuation(x1, y2);
msg = strcat({'Attenuation x12 = '}, num2str(attRT12));
disp(msg{1});
attRT22 = grid_TS.pixelAttenuation(x2, y2);
msg = strcat({'Attenuation x22 = '}, num2str(attRT22));
disp(msg{1});

%==============================
% Filter analysis
%==============================
figure;
plot(grid_TS.tFilter, grid_TS.Filter);
title('Filter')
saveas(gcf, 'Example25_Filter.fig');

%%  % Generate Filter time array
%%  % Set filter characteristics
%%  grid_TS.setCFilter(c0);
%%  grid_TS.setTFilter(tArray_TS);
%%  grid_TS.setFilter(Filter_TS);
%%  % Minimum and maximum sound speed
%%  cMin = 1200;
%%  cMax = 1800;
%%  tMin = grid_TS.source(1).tSignal(2); % time signal increment
%%  aSignal = grid_TS.source(1).aSignal;
%%  lengthASignal = length(aSignal);
%%  tFilter = 0:tMin:grid_TS.tFilter(end)*grid_TS.cFilter/cMax+tMin; % time array for the filter
%%  lengthFilter = length(tFilter);
%%  nFilters = 100; % number of filters;
%%  c = cMin:(cMax - cMin)/nFilters:cMax;
%%  FilterFFTMatrix = zeros(grid_TS.Nx, grid_TS.Ny, lengthFilter + lengthASignal);
%%  cMatrix = repmat(grid_TS.c, [1 1 lengthFilter + lengthASignal]);
%%  for i = 1:length(c)
%%      disp(int2str(i));
%%      % Compute splines for filter
%%      splineFilter = spline(grid_TS.tFilter*grid_TS.cFilter/c(i), grid_TS.Filter, tFilter);
%%      FilterPad = padarray(splineFilter, [0 lengthASignal], 0, 'post');
%%      FilterFFTaux = repmat(permute(fft(FilterPad), [1 3 2]), [grid_TS.Nx grid_TS.Ny 1]);
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
