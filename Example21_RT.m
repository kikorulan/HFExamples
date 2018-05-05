%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

% Measure computational time
tic;
start_time = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 720;           % number of grid points in the x (row) direction
Ny = 1080;           % number of grid points in the y (column) direction
dx = 0.12/Nx;        % grid point spacing in the x direction [m]
dy = 0.18/Ny;        % grid point spacing in the y direction [m]
kernelSize = 33; % Multiple of 3

Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -0.1;
v2 = 0.1;
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

% Build initial pressure
u0 = 1;
M1 = zeros(Nx, floor(Ny/2));
M2 = ones(floor(Nx/3),floor(Ny/4));
M3 = 0*M2;
M4 = zeros(Nx, floor(Ny/4));
U = [M1 [M3; M2; M3] M4];
Rgrid.setUMatrix(U);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays & sources
nSources = 1;

tauMax = sqrt((Nx*dx)^2+(Ny*dy)^2)*c0;
step = min(dx, dy)*c0/3;

% Sources locations
%x1 = [Nx*dx/2; Ny*dy/100]; % Source 1
%x2 = [5*Nx*dx/6; Ny*dy/100]; % Source 2
%x2 = [2*Nx*dx/3; 8*Ny*dy/9]; % Point in the shadow
x3 = [Nx*dx/100; Ny*dy/3]; % Source 3


criticalAngle = 0.7460;
% Create the new sources
nRays = 300;
Rgrid.newSource(x3, criticalAngle-5e-3, criticalAngle+5e-3, nRays, step, tauMax);
%Rgrid.newSource(x3, criticalAngle-5e-1, criticalAngle+5e-1, nRays, step, tauMax);
%Rgrid.newSource(x3, 0, 3*pi/8, nRays, step, tauMax);
%Rgrid.newSource(x3, 0, criticalAngle, nRays, step, tauMax);
%Rgrid.newSource(x2, pi/2, 4*pi/5, nRays, step, tauMax/5);

% Sources
for n = 1:nSources
    disp(strcat('Source ', int2str(n)));
    Rgrid.computeSourceQ(n);
%    Rgrid.revAmplitudeSI(n);
end

nRays = size(Rgrid.source(1).x, 3);

%==================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================================================================================
x = 1:1:Nx;
y = 1:1:Ny;

close all;
run colourMap;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sound speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
flipGrid = Rgrid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_C.fig'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time vs Tay
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Rays
nColours = 8;
colourList = cool(nColours);
hold on;
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(Rgrid.source(1).tau, Rgrid.source(1).phi(:, :, j), 'Color', colourList(colourNum, :));
end
title('Time vs Tau');
xlabel('tau');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_TimeVsTau.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources and their respective Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
axis([0 Nx*dx 0 Ny*dy]);
% Rays
nColours = 8;
colourList = cool(nColours);
hold on;
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(Rgrid.source(1).x(1, :, j), Rgrid.source(1).x(2, :, j), 'Color', colourList(colourNum, :));
end
title('Trajectories');
xlabel('t');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_Ray.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isotime surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCurves = 100;
x = Rgrid.findIsoTime(1, nCurves);
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
% Rays
nColours = 8;
colourList = cool(nColours);
for j = 1:nCurves
     colorNum = floor(nColours*(j-1)/nCurves)+1;
     plot(x(1, :, j), x(2, :, j), 'Color', colourList(colorNum, :));
end
title('Isotime curves');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_Isotime.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude and Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amplitude
figure;
nColours = 8;
colourList = cool(nColours);
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    semilogy(Rgrid.source(1).phi(:, :, j), Rgrid.source(1).amplitude(:, :, j), 'Color', colourList(colourNum, :));
    hold on;
end
grid on;
title('Amplitude');
xlabel('t');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_A.fig');

% Q
figure;
nColours = 8;
colourList = cool(nColours);
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(Rgrid.source(1).phi(:, :, j), Rgrid.source(1).q(:, :, j), 'Color', colourList(colourNum, :));
    hold on;
end
grid on;
title('Q');
xlabel('t');
saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_Q.fig');

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  % Interface loss
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  I = Rgrid.interfaceLoss(1);
%%  figure;
%%  nColours = 8;
%%  colourList = cool(nColours);
%%  for j = 1:nRays
%%      colorNum = mod(floor(nColours*(j-1)/nRays), 6)+1;
%%      semilogy(Rgrid.source(1).phi(:, :, j), I(:, :, j), 'Color', colourMapV(colorNum));
%%      hold on;
%%  end
%%  grid on;
%%  title('Interface Loss');
%%  xlabel('t');
%%  saveas(gcf, '~/Documents/MATLAB/HighFreq/Examples/Ex21_RT_shadow/Example21_Interface.fig');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

