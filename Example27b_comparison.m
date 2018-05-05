% EXAMPLE27B_COMPARISON compares the results between kWave and RT

%load sensorRT;
%load sensorKWAVE;

close all;

% Define colours
darkR = [0.8 0.3 0.3];
darkG = [0.3 0.8 0.3];
darkB = [0.3 0.3 0.8];
darkM = [0.8 0   0.8];

% Plot limits
yMin = -0.7;
yMax = 1.2;
axisY = [0 5e-5 yMin yMax];
yticks = yMin:0.1:yMax;

%==================
% Obtain norms
%==================
% kWave
LinfKW1 = max(sensorKWAVE_angle1.p(2, :));
L2KW1 = sqrt(sum(sensorKWAVE_angle1.p(2, :).^2));
% RT
LinfRT1 = max(sensorRT_angle1(1, :));
L2RT1 = sqrt(sum(sensorRT_angle1(1, :).^2));
% Norms
normKWAVE = LinfKW1;
normRT = LinfRT1;

%====================================================================
% Medium 1
%====================================================================
timeSignalRT = grid_angle1.timeSignal;
timeSignalKWAVE = kgrid.t_array;

%====================
% Superposed plot
%====================
figure;
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(2, :)/normKWAVE, 'Color', darkR, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(3, :)/normKWAVE, 'Color', darkG, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(4, :)/normKWAVE, 'Color', darkB, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(5, :)/normKWAVE, 'Color', darkM, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_angle1(1, :)/normRT, '-r');
plot(timeSignalRT, sensorRT_angle1(2, :)/normRT, '-g');
plot(timeSignalRT, sensorRT_angle1(3, :)/normRT, '-b');
plot(timeSignalRT, sensorRT_angle1(4, :)/normRT, '-m');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW', ...
    'Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensor comparison - Medium 1');
xlabel('t (s)');
ylabel('Amplitude');
grid on;
saveas(gcf, 'Example27b_Sensor_Medium1_A.fig');

%====================
% Subplot
%====================
figure;
subplot(2, 1, 1);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(2, :)/normKWAVE, 'Color', darkR);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(3, :)/normKWAVE, 'Color', darkG);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(4, :)/normKWAVE, 'Color', darkB);
plot(timeSignalKWAVE, sensorKWAVE_angle1.p(5, :)/normKWAVE, 'Color', darkM);
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW');
title('Sensors for kWave - Medium 1');
grid on;
subplot(2, 1, 2);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalRT, sensorRT_angle1(1, :)/normRT, '-r');
plot(timeSignalRT, sensorRT_angle1(2, :)/normRT, '-g');
plot(timeSignalRT, sensorRT_angle1(3, :)/normRT, '-b');
plot(timeSignalRT, sensorRT_angle1(4, :)/normRT, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensors for RT - Medium 1');
grid on;
saveas(gcf, 'Example27b_Sensor_Medium1_B.fig');

%====================================================================
% Medium 2
%====================================================================
timeSignalRT = grid_angle2.timeSignal;
timeSignalKWAVE = kgrid.t_array;

%====================
% Superposed plot
%====================
figure;
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(2, :)/normKWAVE, 'Color', darkR, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(3, :)/normKWAVE, 'Color', darkG, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(4, :)/normKWAVE, 'Color', darkB, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(5, :)/normKWAVE, 'Color', darkM, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_angle2(1, :)/normRT, '-r');
plot(timeSignalRT, sensorRT_angle2(2, :)/normRT, '-g');
plot(timeSignalRT, sensorRT_angle2(3, :)/normRT, '-b');
plot(timeSignalRT, sensorRT_angle2(4, :)/normRT, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW', ...
    'Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensor comparison - Medium 2');
grid on;
saveas(gcf, 'Example27b_Sensor_Medium2_A.fig');

%====================
% Subplot
%====================
figure;
subplot(2, 1, 1);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(2, :)/normKWAVE, 'Color', darkR);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(3, :)/normKWAVE, 'Color', darkG);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(4, :)/normKWAVE, 'Color', darkB);
plot(timeSignalKWAVE, sensorKWAVE_angle2.p(5, :)/normKWAVE, 'Color', darkM);
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW');
title('Sensors for kWave - Medium 2');
grid on;
subplot(2, 1, 2);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalRT, sensorRT_angle2(1, :)/normRT, '-r');
plot(timeSignalRT, sensorRT_angle2(2, :)/normRT, '-g');
plot(timeSignalRT, sensorRT_angle2(3, :)/normRT, '-b');
plot(timeSignalRT, sensorRT_angle2(4, :)/normRT, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensors for RT - Medium 2');
grid on;
saveas(gcf, 'Example27b_Sensor_Medium2_B.fig');

%====================================================================
% Medium 3
%====================================================================
timeSignalRT = grid_angle3.timeSignal;
timeSignalKWAVE = kgrid.t_array;
factorRT = LinfKW1; % max(sensorKWAVE_angle3.p(2, :));
factorKWAVE = L2RT1; % max(sensorRT_angle3(1, :));

%====================
% Superposed plot
%====================
figure;
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(2, :)/normKWAVE, 'Color', darkR, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(3, :)/normKWAVE, 'Color', darkG, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(4, :)/normKWAVE, 'Color', darkB, 'LineWidth', 2);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(5, :)/normKWAVE, 'Color', darkM, 'LineWidth', 2);
plot(timeSignalRT, sensorRT_angle3(1, :)/normRT, '-r');
plot(timeSignalRT, sensorRT_angle3(2, :)/normRT, '-g');
plot(timeSignalRT, sensorRT_angle3(3, :)/normRT, '-b');
plot(timeSignalRT, sensorRT_angle3(4, :)/normRT, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW', ...
    'Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensor comparison - Medium 3');
grid on;
saveas(gcf, 'Example27b_Sensor_Medium3_A.fig');

%====================
% Subplot
%====================
figure;
subplot(2, 1, 1);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(2, :)/normKWAVE, 'Color', darkR);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(3, :)/normKWAVE, 'Color', darkG);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(4, :)/normKWAVE, 'Color', darkB);
plot(timeSignalKWAVE, sensorKWAVE_angle3.p(5, :)/normKWAVE, 'Color', darkM);
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - KW', 'Sensor 2 - KW', 'Sensor 3 - KW', 'Sensor 4 - KW');
title('Sensors for kWave - Medium 3');
grid on;
subplot(2, 1, 2);
hold on;
axis(axisY);
set(gca, 'ytick', yticks);
plot(timeSignalRT, sensorRT_angle3(1, :)/normRT, '-r');
plot(timeSignalRT, sensorRT_angle3(2, :)/normRT, '-g');
plot(timeSignalRT, sensorRT_angle3(3, :)/normRT, '-b');
plot(timeSignalRT, sensorRT_angle3(4, :)/normRT, '-m');
xlabel('t (s)');
ylabel('Amplitude');
legend('Sensor 1 - RT', 'Sensor 2 - RT', 'Sensor 3 - RT', 'Sensor 4 - RT');
title('Sensors for RT - Medium 3');
grid on;
saveas(gcf, 'Example27b_Sensor_Medium3_B.fig');

%====================================================================
% Plot Sound Speeds
%====================================================================
Nx = grid_angle1.Nx;
Ny = grid_angle1.Ny;
dx = grid_angle1.dx;
dy = grid_angle1.dy;
X = 0:dx:(Nx - 1)*dx;
Y = 0:dy:(Ny - 1)*dy;

% Sensor Matrix
sensorMatrix = zeros(Nx, Ny);
sensorMatrix(sensorsRT(1, 1)-2:sensorsRT(1, 1)+2, sensorsRT(2, 1)-2:sensorsRT(2, 1)+2) = inf;
sensorMatrix(sensorsRT(1, 2)-2:sensorsRT(1, 2)+2, sensorsRT(2, 2)-2:sensorsRT(2, 2)+2) = inf;
sensorMatrix(sensorsRT(1, 3)-2:sensorsRT(1, 3)+2, sensorsRT(2, 3)-2:sensorsRT(2, 3)+2) = inf;
sensorMatrix(sensorsRT(1, 4)-2:sensorsRT(1, 4)+2, sensorsRT(2, 4)-2:sensorsRT(2, 4)+2) = inf;
sensorMatrix(sensorsRT(1, 1)-1:sensorsRT(1, 1)+1, sensorsRT(2, 1)-1:sensorsRT(2, 1)+1) = 0;
sensorMatrix(sensorsRT(1, 2)-1:sensorsRT(1, 2)+1, sensorsRT(2, 2)-1:sensorsRT(2, 2)+1) = -300;
sensorMatrix(sensorsRT(1, 3)-1:sensorsRT(1, 3)+1, sensorsRT(2, 3)-1:sensorsRT(2, 3)+1) = -300;
sensorMatrix(sensorsRT(1, 4)-1:sensorsRT(1, 4)+1, sensorsRT(2, 4)-1:sensorsRT(2, 4)+1) = -300;
sensorMatrix(Nx/2-2:Nx/2+2, Ny/10-2:Ny/10+2) = inf;
sensorMatrix(Nx/2-1:Nx/2+1, Ny/10-1:Ny/10+1) = 300;

% Medium 1
figure;
sensorCMatrix = sensorMatrix + grid_angle1.c;
surf(X, Y, sensorCMatrix', 'EdgeColor', 'none');
view(2);
title('Sound Speed & Sensors - Medium 1');
xlabel('x (m)');
ylabel('y (m)');
colorbar();
saveas(gcf, 'Example27b_C_SensorLocations_Medium1.fig');

% Medium 2
figure;
sensorCMatrix = sensorMatrix + grid_angle2.c;
surf(X, Y, sensorCMatrix', 'EdgeColor', 'none');
view(2);
title('Sound Speed & Sensors - Medium 2');
xlabel('x (m)');
ylabel('y (m)');
colorbar();
saveas(gcf, 'Example27b_C_SensorLocations_Medium2.fig');

% Medium 3
figure;
sensorCMatrix = sensorMatrix + grid_angle3.c;
surf(X, Y, sensorCMatrix', 'EdgeColor', 'none');
view(2);
title('Sound Speed - Medium 3');
xlabel('x (m)');
ylabel('y (m)');
colorbar();
saveas(gcf, 'Example27b_C_SensorLocations_Medium3.fig');
