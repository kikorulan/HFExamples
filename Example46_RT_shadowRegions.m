%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex46_caustics;

close all;
%clear all;

load sensor_data.mat;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction
dx = 1e-4;        % grid point spacing in the x direction [m]
dy = 1e-4;        % grid point spacing in the y direction [m]
gridR = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = -.2;
v2 = 0.1;
kernelSize = 20;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
c = M1;
%c = addCircle(c, kernelSize + floor(Nx/2), kernelSize + floor(Ny/2), floor(Nx/6), c0*v1);
%c = addCircle(c, kernelSize + floor(Nx-1), kernelSize + floor(Ny/2), floor(Nx/2), c0*v1);
c(floor(5*Nx/9):end, :) = c(floor(5*Nx/9):end, :) + c0*v1;
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.density = 1;

c = medium.sound_speed;
gridR.setCMatrix(c);

%========================================
% Impulse Response
%========================================
% Set time
dt = 2e-8;
%dt = min(gridR.dx, gridR.dy)/c0/2;
tMax = 2.5e-5;
gridR.setTime(dt, tMax);
% Compute impulse response
gridR.impulse_additive('IV');
save gridRT_impulse_shadow.mat gridR c0 dt tMax;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse_shadow.mat

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 800;% 800
nSources = 1;%256

% Parametrisation
tStep = dt;

% Sources locations
clear x;
x{1} = cat(3, (gridR.Nx-1)/2*gridR.dx, 0);

% Sources
source(1) = gridR.newSource(x{1}, 0, pi, nRays, tStep, tMax);
% Set initial pressure
gridR.setUMatrix(source_low.p0);
gridR.computeHamil(source(1));

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex46_caustics;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

%load gridRRT.mat;
%load sensor_data.mat;
position = [700 500 320 600];
positionYBar = [700 700 375 600];
positionBar = [700 700 375 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Sound speed
%==============================
gridR.plot_soundSpeed();
set(gcf, 'pos', positionYBar);
colorbar();
saveas(gcf, 'Example46_C_shadow', 'png'); 
saveas(gcf, 'Example46_C_shadow.fig'); 

%==============================
% Traveled distance
%==============================
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      plot(gridR.tForward, real(source(1).xTD(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Traveled distance');
%%  xlabel('t (s)');

%==============================
% Sources and their respective Rays
%==============================
%%  for n = 1:nSources
%%     source(1).plot_rays(gridR, 10);
%%  end

%==============================
% Isotime curve
%==============================
%%  nCurves = 20;
%%  xIso = gridR.findIsoTime(source(1), 20, 1e-4, 1.2e-4);
%%  % Figure
%%  figure;
%%  flipGrid = 1./gridR.c';
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  hold on;
%%  surf(gridR.xAxis, gridR.yAxis, flipGrid, 'EdgeColor', 'none');
%%  hold on;
%%  for n = 1:nCurves
%%      plot3(xIso(n, :, 1), xIso(n, :, 2), repmat(1501, [1 size(xIso, 2)]), '-m');
%%  end
%%  for j = 1:20:nRays
%%      plot3(source(2).x(j, :, 1), source(2).x(j, :, 2), repmat(1501, [1 source(2).nPoints]), '-g');
%%  end
%%  view(2);
%%  legend('Isocurves for times tMin = 1e-4s, tMax = 1.2e-4s');

%==============================
% Sound speed + rays
%==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = nRays;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:nRays
%%      colourNum = floor(nColours*(j-1)/nRays) + 1;
%%      plot3(source(1).x(j, :, 1), source(1).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  %surf(gridR.xAxis, gridR.yAxis, 1./gridR.c', 'EdgeColor', 'none');
%%  imagesc(gridR.xAxis, gridR.yAxis, 1./gridR.c');
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', position);
%%  %saveas(gcf, 'Example46_rays_soundSpeed.fig');

%==============================
% Ray trajectories
%==============================
h = figure;
hold on;
axis([0 1e3*gridR.xAxis(end) 0 1e3*gridR.yAxis(end)]);
nColours = nRays;
colourList = cool(nColours);
[nRays nSteps dim] = size(source(1).x);
zVec = 1600*ones(1, nSteps);
for j = 1:5:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot3(1e3*source(1).x(j, :, 1), 1e3*source(1).x(j, :, 2), zVec, 'Color', colourList(colourNum, :));
end
xlabel('x [mm]');
ylabel('y [mm]');
box on;
set(gcf, 'pos', position);
saveas(gcf, 'Example46_rays_shadow', 'png');
saveas(gcf, 'Example46_rays_shadow.fig');


%==============================
% Amplitude
%==============================
%%  % Imagesc
%%  figure;
%%  surf(real(log(real(source(1).amplitude))), 'EdgeColor', 'none');
%%  view(2);
%%  xlabel('time step');
%%  ylabel('ray number');
%%  title('Amplitudes (log scale)');
%%  colorbar();
%%  hold on;
%%  %plot(1:1000, 350*ones(1, 1000), 'Color', 'r');
%%  %plot(1:1000, 450*ones(1, 1000), 'Color', 'r');
%%  %%  saveas(gcf, 'Example46_surf_amplitudeDin_log', 'png');
%%  %%  
%%  % Log plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      semilogy(gridR.tForward, real(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Amplitudes for rays with caustic (Gaussian Beam)');
%%  xlabel('t (s)');
%%  saveas(gcf, 'Example46_amplitudeGB_log', 'epsc');
%%  
%%  % Log plot
%%  figure;
%%  colours = winter(100);
%%  for n = 1:100
%%      semilogy(gridR.tForward, imag(source(1).amplitude(n+350, :)), 'Color', colours(n, :));
%%      hold on;
%%  end;
%%  grid on;
%%  box on;
%%  title('Amplitudes for rays with caustic - imag');
%%  xlabel('t (s)');
%%  %saveas(gcf, 'Example46_amplitudeODE_log_caustic.fig');



%==============================
% Sound speed + ray selection
%==============================
%%  h = figure;
%%  hold on;
%%  axis([0 gridR.xAxis(end) 0 gridR.yAxis(end)]);
%%  nColours = 320;
%%  colourList = summer(nColours);
%%  [nRays nSteps dim] = size(source(1).x);
%%  zVec = 1600*ones(1, nSteps);
%%  for j = 1:5:nColours
%%      colourNum = j;
%%      plot3(source(1).x(j + 400-nColours/2, :, 1), source(1).x(j + 400-nColours/2, :, 2), zVec, 'Color', colourList(colourNum, :));
%%  end
%%  c = parula;
%%  c = flipud(c);
%%  colormap(c);
%%  imagesc(gridR.xAxis, gridR.yAxis, gridR.c');
%%  colorbar();
%%  xlabel('x (m)');
%%  ylabel('y (m)');
%%  title('Ray Trajectories - focus');
%%  set(gcf, 'pos', positionYBar);
%%  saveas(gcf, 'Example46_rays_soundSpeed_selection', 'png');
%%  saveas(gcf, 'Example46_rays_soundSpeed_selection.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;


