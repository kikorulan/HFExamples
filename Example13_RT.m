%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for RT - Analyse impact of deltaX and deltaP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear grid A B;

% Run color map
run colourMap;
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

% Set deltaX and deltaP
grid.setDeltaX(grid.dx/1e3);
grid.setDeltaP(0.3);
% Build domain
c0 = 1500;
% Auxiliary matrices
c = c0*ones(Nx, Ny);
grid.setCMatrix(c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays
nRays = 1;
nSources = 3;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
x1 = [Nx*dx/2; Ny*dy/100]; % Source 1

% Sources
grid.setDeltaP(0.3);
% Study Delta X variations 
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    % Create the new sources
    grid.setDeltaX(grid.dx*10^(-n));
    grid.newSource(x1, pi/2, pi/2, nRays, step, tauMax);
    grid.computeSource(n);
    A{n} = grid.computeAmplitudeODE(n);
    B{n} = grid.computeAmplitudeProximal(n);
%    grid.computeReverseAmplitude(n);
end

% Sources
grid.setDeltaX(grid.dx/1e2);
% Study Delta P variations
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    % Create the new sources
    grid.setDeltaP(10^(-n/2));
    grid.newSource(x1, pi/2, pi/2, nRays, step, tauMax);
    grid.computeSource(n+nSources);
    A{n+nSources} = grid.computeAmplitudeODE(n+nSources);
    B{n+nSources} = grid.computeAmplitudeProximal(n+nSources);
%    grid.computeReverseAmplitude(n);
end


%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude - Delta X
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Sources
for n = 1:nSources
    phi = grid.source(n).phi(:, :, 1);
    loglog(phi(1:length(A{n})), A{n}, 'Color', colourMapRed(n), 'LineWidth', 4);
    hold on;
    loglog(phi(1:length(B{n})), B{n}, 'Color', colourMapBlue(n), 'LineWidth', 0.01, ...
             'Marker', 'o', 'MarkerSize', 1);
end
% Square root
loglog(phi(phi>0), 1./sqrt(phi(phi>0)), '-xg', 'MarkerSize', 3);
legend('ODE - deltaX = grid.dx*1e-1', 'PrR - deltaX = grid.dx*1e-1', ...
       'ODE - deltaX = grid.dx*1e-2', 'PrR - deltaX = grid.dx*1e-2', ...
       'ODE - deltaX = grid.dx*1e-3', 'PrR - deltaX = grid.dx*1e-3', ...
       '1/sqrt(t)');
title('Delta X for ODE and Proximal Ray');
xlabel('t');
ylabel('Amplitude');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex13_RT_deltaX/Example13_Amplitude_DeltaX.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude - Delta P
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Sources
phi = grid.source(nSources+1).phi(:, :, 1);
loglog(phi(1:length(A{nSources + 1})), A{nSources + 1}, 'Color', colourMapRed(n), 'LineWidth', 4);
hold on;
for n = 1:nSources
    phi = grid.source(n+nSources).phi(:, :, 1);
    loglog(phi(1:length(B{n+nSources})), B{n+nSources}, 'Color', colourMapBlue(n), 'LineWidth', 0.01, ...
             'Marker', 'o', 'MarkerSize', 1);
end
% Square root
loglog(phi(phi>0), 1./sqrt(phi(phi>0)), '-xg', 'MarkerSize', 3);
legend('ODE', 'PrR - deltaP = 1e-1', ...
       'PrR - deltaP = 1e-2', ...
       'PrR - deltaP = 1e-3', ...
       '1/sqrt(t)');
title('Delta P for ODE and Proximal Ray');
xlabel('t');
ylabel('Amplitude');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex13_RT_deltaX/Example13_Amplitude_DeltaP.fig'); 

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

