%================================================================================
%% Example for gridRT class
%================================================================================
cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex21_shadow;

close all;
%clear all;

%load('workspace');
% Measure computational time
tic;
start_time = clock;

run colourMap;

%========================================
% Grid definition
%========================================
Nx = 720;            % number of grid points in the x (row) direction
Ny = 1080;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
kernelSize = 33; % Multiple of 3

Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.2;
v2 = 0.2;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, floor(dimY/3));
M2 = v1*c0*ones(floor(dimX/2), floor(dimY/3));
M3 = v2*c0*ones(floor(dimX/2), floor(dimY/3));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [M1...
     [M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
Rgrid.setCMatrix(cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize));

%========================================
% Ray Shooting
%========================================

% Number of rays & sources
nRays = 100;
nRaysSub = 30;
nSources = 4;

% Tau parameter
tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
tauStep = min(dx, dy)*c0/3;
% T parameter
tMax = sqrt((Nx*dx)^2+(Ny*dy)^2)/c0;
tStep = min(dx, dy)/c0/3;

% Sources locations
x2 = cat(3, 5*Nx*dx/6, Ny*dy/100); % Source 2
% Create the new sources
%Rgrid.newSource(x2, pi/2, 5*pi/6, nRays, tauStep, tauMax);
% Subdivide rays tau-Hamiltonian
Rgrid.newSource(x2, pi/2, 13*pi/24, nRaysSub, tauStep, tauMax);
Rgrid.newSource(x2, 13*pi/24, 43*pi/72, nRaysSub, tauStep, tauMax);
Rgrid.newSource(x2, 43*pi/72, 5*pi/6, nRaysSub, tauStep, tauMax);
% Subdivide rays t-Hamiltonian
Rgrid.newSource(x2, pi/2, 13*pi/24, nRaysSub, tStep, tMax);
Rgrid.newSource(x2, 13*pi/24, 43*pi/72, nRaysSub, tStep, tMax);
Rgrid.newSource(x2, 43*pi/72, 5*pi/6, nRaysSub, tStep, tMax);


% Sources
t1 = clock;
for n = 1:3
    disp(strcat('Source ', int2str(n)));
    Rgrid.computeTrajec(n);
    Rgrid.amplitudeProximal(n);
end

for n = 4:6
    disp(strcat('Source ', int2str(n)));
    Rgrid.computeHamil(n);
end


t2 = clock;
disp(['  Computation time RT ' num2str(etime(t2, t1))]);

save gridRT.mat Rgrid;

%==================================================================================
% Plot results
%==================================================================================
close all;
load gridRT.mat

cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex21_shadow;

x = 0:Rgrid.dx:(Rgrid.Nx-1)*Rgrid.dx;
y = 0:Rgrid.dy:(Rgrid.Ny-1)*Rgrid.dy;
  
%%  %==============================
%%  % Sound speed
%%  %==============================
%%  figure;
%%  flipGrid = Rgrid.c';
%%  surf(x, y, flipGrid, 'EdgeColor', 'none');
%%  view(2);
%%  colorbar();
%%  title('Sound Speed');
%%  saveas(gcf, 'Example21b_C.fig'); 
%%  
%%  %==============================
%%  % Sources and their respective Rays
%%  %==============================
%%  Rgrid.plotRays([1 2 3], nRays);
%%  title('Trajectories');
%%  saveas(gcf, 'Example21b_RaySub_tauH.fig');
%%  
%%  Rgrid.plotRays([4 5 6], nRays);
%%  title('Trajectories');
%%  saveas(gcf, 'Example21b_RaySub_tH.fig');

%==============================
% Amplitude and Q
%==============================
% Amplitude
for n = 1:3
    Rgrid.plotAmplitude(n, 100);
    title('Amplitude');
    %saveas(gcf, strcat('Example21b_RayA', int2str(n), '_tH.fig'));
end

for n = 4:6
    Rgrid.plotAmplitude(n, 100);
    title('Amplitude');
    %saveas(gcf, strcat('Example21b_RayA', int2str(n), '_tauH.fig'));
end

% Q and rays
for n = 1:3
    figure;
    colourList = cool(nRaysSub);
    % Subplot 1
    subplot(1, 2, 1);
    axis([0 Rgrid.Nx*Rgrid.dx 0 Rgrid.Ny*Rgrid.dx]);
    hold on;
    for j = 1:nRaysSub
        plot(Rgrid.source(n).x(j, :, 1), Rgrid.source(n).x(j, :, 2), 'Color', colourList(j, :));
    end
    title('Trajectories');
    xlabel('x');
    ylabel('y');
    grid on;
    % Subplot 2
    subplot(1, 2, 2);
    hold on;
    for j = 1:nRaysSub
        plot(Rgrid.source(n).phi(j, :), Rgrid.source(n).q(j, :), 'Color', colourList(j, :));
    end
    title('Q - tauH');
    xlabel('t (s)');
    grid on;
    %saveas(gcf, strcat('Example21b_RayQ', int2str(n), '_tauH.fig'));
end

for n = 4:6
    figure;
    colourList = cool(nRaysSub);
    % Subplot 1
    subplot(1, 2, 1);
    axis([0 Rgrid.Nx*Rgrid.dx 0 Rgrid.Ny*Rgrid.dx]);
    hold on;
    for j = 1:nRaysSub
        plot(Rgrid.source(n).x(j, :, 1), Rgrid.source(n).x(j, :, 2), 'Color', colourList(j, :));
    end
    title('Trajectories');
    xlabel('x');
    ylabel('y');
    grid on;
    % Subplot 2
    subplot(1, 2, 2);
    hold on;
    for j = 1:nRaysSub
        plot(Rgrid.source(n).phi(j, :), Rgrid.source(n).q(j, :), 'Color', colourList(j, :));
    end
    title('Q - tH');
    xlabel('t (s)');
    grid on;
    %saveas(gcf, strcat('Example21b_RayQ', int2str(n), '_tH.fig'));
end
 
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

