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

% Set deltaX and deltaP
grid.setDeltaX(grid.dx/1e3);
grid.setDeltaP(0.3);
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
nRays = 1;
nSources = 1;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
x1 = [Nx*dx/2; Ny*dy/100]; % Source 1
x2 = [5*Nx*dx/6; Ny*dy/100]; % Source 2
x3 = [Nx*dx/100; Ny*dy/3]; % Source 3

% Create the new sources
grid.newSource(x1, pi/3, 2*pi/3, nRays, step, tauMax);
grid.newSource(x2, pi/2, 5*pi/6, nRays, step, tauMax);
grid.newSource(x3, 0, 3*pi/8, nRays, step, tauMax);

% Sources
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    grid.computeSource(n);
    A = grid.computeAmplitude_ODE(n);
    B = grid.computeAmplitude_ODEct(n);
    grid.computeAmplitude(n,1);
%    grid.computeReverseAmplitude(n);
end

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;


%%%%
% Sources and their respective Rays
%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
% Sources
for n = 1:nSources
    % Rays
    for j = 1:nRays
        plot(grid.source(n).x(1, :, j), grid.source(n).x(2, :, j), 'Color', colourMap(n));
    end
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex12_RT_ODE/Example12_Ray.fig');


%%%%
% Phase
%%%%
figure;
surf(x, y, grid.phase, 'EdgeColor', 'none');
view(2);


% Amplitude
figure;
hold on;
% Sources
for n = 1:nSources
    % Rays
    for j = 1:nRays
        ray = grid.source(n).amplitude(:, :, j);
        phi = grid.source(n).phi(:, :, j);
        plot(phi(1:length(A)), log(A), 'Color', colourMap(n));
        plot(phi(1:length(B)), log(B), 'Color', colourMap(n+1));
        plot(phi(1:length(ray)), log(ray), 'Color', colourMap(n+2));
    end
end
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex12_RT_ODE/Example12_Amplitude.fig'); 
legend('ODE', 'ODE constant c', 'Ray');


end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

