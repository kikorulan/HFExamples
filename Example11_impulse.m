%================================================================================
% Example for gridRT class
%================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex11_source;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex11_source;

close all;
%clear all;

%load gridRT_ray.mat;
%load kgrid_data.mat;

% Measure computational time
tic;
start_time = clock;

run colourMap;
%========================================
% Impulse Response
%========================================
grid.impulseResponse2D('GF');

%========================================
% Forward Signal
%========================================

tSpline = 0:grid.dt/10:grid.tForward(end)+grid.dt/10;
lForward = length(tSpline);
for n = 1:3
    % Generate Spline
    aSpline = spline(grid.source(n).tBeam, grid.source(n).aBeam, tSpline);
    % Convolve Signal and Filter
    signalConv = conv(aSpline, grid.Filter);
    % Adjust delay
    delayFilter = find(grid.Filter == max(grid.Filter)) - 1;
    aForwardDelay = signalConv(1+delayFilter:end);
    % Adjust total length
    if (length(aForwardDelay) < lForward)
        aForward = padarray(aForwardDelay, [0 lForward-length(aForwardDelay)], 0, 'post');
    else 
        aForward = aForwardDelay(1:lForward);
    end
    % Save Signal
    grid.setForwardSignal(n, aForward);
end


%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex11_source;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex11_source;

%load gridRT.mat;
load kgrid_data.mat;

%==============================
% Non-Filtered Signal
%==============================
normBeam = max(grid.source(1).aBeam);
% Time Signal Amplitude
figure;
hold on;
for n = 1:nSources
    plot(grid.source(n).tBeam, grid.source(n).aBeam/normBeam, 'Color', colourMapV(n));  
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3')

%==============================
% Filtered Signal
%==============================
% Normalisation factor
normRT = max(grid.source(1).aForward);
normKWave = max(sensor_data.p(1, :));

% Figure
figure;
hold on;
%axis([0 1.2e-4 -0.5 0.5]);
for n = 1:nSources
    % Plot and compare
    plot(tSpline, grid.source(n).aForward/normRT, 'Color', colourMapV(n));
    plot(kgrid.t_array, sensor_data.p(n, :)/normKWave, 'Color', colourMapV(n+nSources), 'LineWidth', 2);
end

legend('Sensor 1 RT', 'Sensor 1 kWave', ...
    'Sensor 2 RT', 'Sensor 2 kWave', 'Sensor 3 RT', 'Sensor 3 kWave');

% Filter
figure;
plot(grid.tFilter, grid.Filter);
title('Filter');
%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/;


