% EXAMPLE27B_COMPARISON compares the results between kWave and RT
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex64_GB;

%load sensorRT;

close all;

% Define colours
darkR = [0.8 0.3 0.3];
darkG = [0.3 0.8 0.3];
darkB = [0.3 0.3 0.8];
darkM = [0.8 0   0.8];

% Plot limits
yMin = -1;
yMax = 1.2;
axisY = [0 5e-5 yMin yMax];
yticks = yMin:0.1:yMax;

%==================
% Obtain norms
%==================
% kWave
LinfKW1 = max(sensorKWAVE.p(1, :));
L2KW1 = sqrt(sum(sensorKWAVE.p(1, :).^2));
% RT
LinfRT1 = max(sensorRT_s1.aForward);
L2RT1 = sqrt(sum(sensorRT_s1.aForward.^2));
% GB
LinfGB1 = max(sensorGB_s1.aForward);
L2GB1 = sqrt(sum(sensorGB_s1.aForward.^2));
% Norms
normKWAVE = LinfKW1;
normRT = LinfRT1;
normGB = LinfGB1;

%====================================================================
% RT vs GB
%====================================================================
timeSignalRT = grid.tForward;
timeSignalGB = grid.tForward;

%====================
% Superposed plot
%====================
figure;
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalRT, sensorRT_s1.aForward/normRT, 'Color', darkR, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s2.aForward/normRT, 'Color', darkG, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s3.aForward/normRT, 'Color', darkB, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s4.aForward/normRT, 'Color', darkM, 'LineWidth', 2);
plot(timeSignalGB, sensorGB_s1.aForward/normGB, '-r');
plot(timeSignalGB, sensorGB_s2.aForward/normGB, '-g');
plot(timeSignalGB, sensorGB_s3.aForward/normGB, '-b');
plot(timeSignalGB, sensorGB_s4.aForward/normGB, '-m');
legend('Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT', ...
    'Sensor 1 - GB', 'Sensor 2 - GB', 'Sensor 3 - GB', 'Sensor 4 - GB');
title('Sensor comparison');
xlabel('t (s)');
ylabel('Amplitude');
grid on;
saveas(gcf, 'Example_Sensor_RTvsGB.fig');

%====================
% Subplot
%====================
figure;
subplot(2, 1, 1);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalRT, sensorRT_s1.aForward/normRT, 'Color', darkR, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s2.aForward/normRT, 'Color', darkG, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s3.aForward/normRT, 'Color', darkB, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s4.aForward/normRT, 'Color', darkM, 'LineWidth', 2);

xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensors for RT');
grid on;
subplot(2, 1, 2);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalGB, sensorGB_s1.aForward/normGB, '-r');
plot(timeSignalGB, sensorGB_s2.aForward/normGB, '-g');
plot(timeSignalGB, sensorGB_s3.aForward/normGB, '-b');
plot(timeSignalGB, sensorGB_s4.aForward/normGB, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - GB', 'Sensor 2 - GB', 'Sensor 3 - GB', 'Sensor 4 - GB');
title('Sensors for GB');
grid on;
saveas(gcf, 'Example_Sensor_RTvsGB_subplot.fig');

%====================
% Error
%====================
error1 = sensorRT_s1.aForward/normRT - sensorGB_s1.aForward/normGB;
error2 = sensorRT_s2.aForward/normRT - sensorGB_s2.aForward/normGB;
error3 = sensorRT_s3.aForward/normRT - sensorGB_s3.aForward/normGB;
error4 = sensorRT_s4.aForward/normRT - sensorGB_s4.aForward/normGB;

figure;
hold on;
%set(gca, 'ytick', yticks);
plot(timeSignalGB, error1, '-r');
plot(timeSignalGB, error2, '-g');
plot(timeSignalGB, error3, '-b');
plot(timeSignalGB, error4, '-m');
legend('Sensor 1 - error', 'Sensor 2 - error', 'Sensor 3 - error', 'Sensor 4 - error');
title('Sensor comparison');
xlabel('t (s)');
ylabel('Amplitude');
box on;
saveas(gcf, 'Example64_error_RTvsGB.fig');

%==================================================================================
% Amplitude
%==================================================================================
%%  % Sensor 1
%%  figure;
%%  surf(amplitude{1}, 'EdgeColor', 'none');
%%  view(2);
%%  figure;
%%  surf(amplitudeGB{1}, 'EdgeColor', 'none');
%%  view(2);
%%  
%%  % Sensor 2
%%  figure;
%%  surf(amplitude{2}, 'EdgeColor', 'none');
%%  view(2);
%%  figure;
%%  surf(amplitudeGB{2}, 'EdgeColor', 'none');
%%  view(2);
%%  
%%  % Sensor 3
%%  figure;
%%  surf(amplitude{3}, 'EdgeColor', 'none');
%%  view(2);
%%  figure;
%%  surf(amplitudeGB{3}, 'EdgeColor', 'none');
%%  view(2);
%%  
%%  % Sensor 4
%%  figure;
%%  surf(amplitude{4}, 'EdgeColor', 'none');
%%  view(2);
%%  figure;
%%  surf(amplitudeGB{4}, 'EdgeColor', 'none');
%%  view(2);

%==================================================================================================
% Plot Sound Speed
%==================================================================================================
positionYBar = [700 700 390 870];
set(0,'DefaultFigurePaperPositionMode','auto');
axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

Nx = grid.Nx;
Ny = grid.Ny;
dx = grid.dx;
dy = grid.dy;
X = 0:dx:(Nx - 1)*dx;
Y = 0:dy:(Ny - 1)*dy;

% Define the sensor
sensor1 = [Nx/2; 2*Ny/10];
sensor2 = [Nx/2; Ny-3*Ny/10];
sensor3 = [Nx/2; Ny-2*Ny/10];
sensor4 = [Nx/2; Ny-1*Ny/10];
sensorsRT = [sensor1 sensor2 sensor3 sensor4];

% Sensor Matrix
sensorMatrix = zeros(Nx, Ny);
sensorMatrix(sensorsRT(1, 1)-2:sensorsRT(1, 1)+2, sensorsRT(2, 1)-2:sensorsRT(2, 1)+2) = inf;
sensorMatrix(sensorsRT(1, 2)-2:sensorsRT(1, 2)+2, sensorsRT(2, 2)-2:sensorsRT(2, 2)+2) = inf;
sensorMatrix(sensorsRT(1, 3)-2:sensorsRT(1, 3)+2, sensorsRT(2, 3)-2:sensorsRT(2, 3)+2) = inf;
sensorMatrix(sensorsRT(1, 4)-2:sensorsRT(1, 4)+2, sensorsRT(2, 4)-2:sensorsRT(2, 4)+2) = inf;
sensorMatrix(sensorsRT(1, 1)-1:sensorsRT(1, 1)+1, sensorsRT(2, 1)-1:sensorsRT(2, 1)+1) = 0;
sensorMatrix(sensorsRT(1, 2)-1:sensorsRT(1, 2)+1, sensorsRT(2, 2)-1:sensorsRT(2, 2)+1) = -150;
sensorMatrix(sensorsRT(1, 3)-1:sensorsRT(1, 3)+1, sensorsRT(2, 3)-1:sensorsRT(2, 3)+1) = -150;
sensorMatrix(sensorsRT(1, 4)-1:sensorsRT(1, 4)+1, sensorsRT(2, 4)-1:sensorsRT(2, 4)+1) = -150;

% Medium 
figure;
sensorCMatrix = sensorMatrix + grid.c;
surf(X, Y, sensorCMatrix', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
title('Sound Speed & Sensors');
xlabel('x (m)');
ylabel('y (m)');
colorbar();
box on;
set(gcf, 'pos', positionYBar);
saveas(gcf, 'Example64_Sensors.fig');

%==================================================================================================
% Initial pressure
%==================================================================================================
figure;
surf(grid.xAxis, grid.yAxis, grid.u', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
title('Initial Pressure');
xlabel('x (m)');
ylabel('y (m)');
colorbar();
box on;
set(gcf, 'pos', positionYBar);
saveas(gcf, 'Example64_InitialPressure.fig');

%====================================================================
% Medium 1
%====================================================================
timeSignalRT = grid.tForward;
timeSignalKWAVE = kgrid.t_array;

%====================
% Superposed plot
%====================
figure;
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE.p(1, :)/normKWAVE, 'Color', darkR, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE.p(2, :)/normKWAVE, 'Color', darkG, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE.p(3, :)/normKWAVE, 'Color', darkB, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE.p(4, :)/normKWAVE, 'Color', darkM, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_s1.aForward/normRT, '-r');
plot(timeSignalRT, sensorRT_s2.aForward/normRT, '-g');
plot(timeSignalRT, sensorRT_s3.aForward/normRT, '-b');
plot(timeSignalRT, sensorRT_s4.aForward/normRT, '-m');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW', ...
    'Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensor comparison - Medium 1');
xlabel('t (s)');
ylabel('Amplitude');
%grid on;
%saveas(gcf, 'Example27b_Sensor_Medium1_A.fig');

%====================
% Subplot
%====================
figure;
subplot(2, 1, 1);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE.p(1, :)/normKWAVE, 'Color', darkR);
plot(timeSignalKWAVE, sensorKWAVE.p(2, :)/normKWAVE, 'Color', darkG);
plot(timeSignalKWAVE, sensorKWAVE.p(3, :)/normKWAVE, 'Color', darkB);
plot(timeSignalKWAVE, sensorKWAVE.p(4, :)/normKWAVE, 'Color', darkM);
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW');
title('Sensors for kWave - Medium 1');
%grid on;
subplot(2, 1, 2);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalRT, sensorRT_s1.aForward/normRT, '-r');
plot(timeSignalRT, sensorRT_s2.aForward/normRT, '-g');
plot(timeSignalRT, sensorRT_s3.aForward/normRT, '-b');
plot(timeSignalRT, sensorRT_s4.aForward/normRT, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensors for RT - Medium 1');
%grid on;
%saveas(gcf, 'Example27b_Sensor_Medium1_B.fig');
%====================
% Error
%====================
error1 = sensorRT_s1.aForward/normRT - sensorGB_s1.aForward/normGB;
error2 = sensorRT_s2.aForward/normRT - sensorGB_s2.aForward/normGB;
error3 = sensorRT_s3.aForward/normRT - sensorGB_s3.aForward/normGB;
error4 = sensorRT_s4.aForward/normRT - sensorGB_s4.aForward/normGB;

figure;
hold on;
%set(gca, 'ytick', yticks);
plot(timeSignalGB, error1, '-r');
plot(timeSignalGB, error2, '-g');
plot(timeSignalGB, error3, '-b');
plot(timeSignalGB, error4, '-m');
legend('Sensor 1 - error', 'Sensor 2 - error', 'Sensor 3 - error', 'Sensor 4 - error');
title('Sensor comparison');
xlabel('t (s)');
ylabel('Amplitude');
box on;
saveas(gcf, 'Example64_error_RTvsGB.fig');

