%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for RT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear grid;

% Run color map
run colorMap;
% Measure computational time
tic;
start_time = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 120;
dx = 1e-3;
Ny = 120;
dy = 1e-3;
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
% Auxiliary matrices
c = c0*ones(Nx, Ny);
grid.setCMatrix(c);

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];
grid.setUMatrix(U);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays
nRays = 100;
nSources = 3;

tauMax = 2*(Nx*dx+Ny*dy)*c0;
step = min(dx, dy)*c0/12;

% Sources locations
x1 = [Nx*dx/2; Ny*dy/100]; % Source 1
x2 = [5*Nx*dx/6; Ny*dy/100]; % Source 2
x3 = [Nx*dx/100; Ny*dy/3]; % Source 3

% Create the new sources
grid.newSource(x1, step, tauMax);
grid.newSource(x2, step, tauMax);
grid.newSource(x3, step, tauMax);

% Sources
thetaMin = [pi/3 pi/2 0]; % Shooting angles
thetaMax = [2*pi/3 5*pi/6 3*pi/8]; % Shooting angles
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    % Rays
    for j = 1:nRays
        theta = thetaMin(n) + (j-1)*(thetaMax(n)-thetaMin(n))/nRays;
        p0 = [cos(theta); sin(theta)]/c0;
        grid.computeRay(n, p0);
        % Reverse Amplitude
        grid.computeReverseAmplitude(n, j);
    end
end
%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%
% N
%%%%
figure;
flipGrid = grid.u'; % transpose
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_U.fig');


%%%%
% Sources and their respective Rays
%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
% Rays
for j = 1:nRays
    % Sources
    for n = 1:nSources
        plot(grid.getX(n, j, 'xCoord'), grid.getX(n, j, 'yCoord'), 'Color', colourMap(n));
        plot(grid.getXPhi(n, j, 'xCoord'), grid.getXPhi(n, j, 'yCoord'), 'Color', colourMap(n));
    end
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_Ray.fig');

% Amplitude
figure;
hold on;
% Sources
for n = 1:nSources
    % Rays
    for j = 1:nRays
        ray = grid.getRevA(n, j);
        phi = grid.getPhi(n, j);
        plot(phi(1:length(ray)), ray, 'Color', colourMap(n));  
    end
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_RevAmplitude.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Non-Filtered Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Signal Amplitude
figure;
hold on;
for n = 1:nSources
    grid.timeSignal(n);
    plot(grid.source(n).tSignal, grid.source(n).aSignal, 'Color', colourMap(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_aSignal.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Smooth Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lsmooth = 19; % only odd values
smoother = 1/Lsmooth*ones(1, Lsmooth);
figure;
hold on;
for n = 1:nSources
    smoothSignal{n} = conv(grid.source(n).aSignal, smoother);
    smoothSignal{n} = smoothSignal{n}(1+floor(Lsmooth/2):end-floor(Lsmooth/2));
    plot(grid.source(n).tSignal, smoothSignal{n}, 'Color', colourMap(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_aSignal_smooth.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Filtered Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Spline
inc = 7.7e-8;
tSpline = 0:inc:1e-4;
% Generate filter
t0 = 7.78e-8;
c = 1.5e3;
t = 100*t0:inc:1e4*t0;
G = 1./((t.*t)-t0*t0).^(0.5)/c;
delta = [1 -1];
FilterSinc = conv(delta, G);
FilterSinc = FilterSinc(1:700);
FilterSinc(1) = -sum(FilterSinc(2:end));
% Figure
figure;
hold on;
for n = 1:nSources
    timeL = length(tSpline);
    %sConv = conv(grid.source(n).aSignal, Filter);
    sSpline{n} = spline(grid.source(n).tSignal, grid.source(n).aSignal, tSpline);
    n
    sConv = conv(sSpline{n}, FilterSinc);
    plot(tSpline, sConv(1:timeL)*inc*4e3, 'Color', colourMap(n));
    plot(kgrid.t_array, sensor_data.p(n, :), 'Color', colourMap(n+nSources));
end
legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_aSignalConv.fig');


end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);



end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Deconvolve Signal
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Time Signal Amplitude
% figure;
% hold on;
% for n = 1:nSources
%     timeL = length(grid.source(n).tSignal);
%     sensorDataSpline{n} = spline(kgrid.t_array, sensor_data.p(n, :), grid.source(n).tSignal);
%     fftSensor{n} = fft(sensorDataSpline{n});
%     fftSignal{n} = fft(smoothSignal{n});
%     h = ifft(fftSensor{n}./fftSignal{n});
%     plot(grid.source(n).tSignal, h, 'Color', colourMap(n));
% end
% legend('h1', 'h2', 'h3');
% saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex10_RT_kWave/Example10_aSignalConv.fig');
 
