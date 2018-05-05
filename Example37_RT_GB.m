%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex37_gaussBeams;

close all;
%clear all;

run colourMap;

%========================================
% Grid definition
%========================================
Nx = 128;         % number of grid points in the x (row) direction
Ny = 256;         % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
grid = gridRT(Nx, dx, Ny, dy);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1500 ; %1500
v1 = 0; %0.1
v2 = 0;
kernelSize = 30;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
M2 = v1*c0*ones(dimX, floor(dimY/2));
M3 = v2*c0*ones(dimX, floor(dimY/2));
c = M1 + [M2 M3]; 
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
% Build domain
grid.setCMatrix(sound_speed);

%%  %========================================
%%  % Impulse Response
%%  %========================================
%%  % Set time
    cMax = max(grid.c(:));
    dt = min(grid.dx, grid.dy)/cMax/2;
%%  dt = 2.5e-7;
%%  tMax = 2e-4;
%%  grid.setTime(dt, tMax);
%%  % Compute impulse response
%%  grid.impulseResponse2D('IV');
%%  
%%  save gridRT_impulse.mat grid;

%========================================
% Ray Shooting
%========================================
load gridRT_impulse;
% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 100;% 800

% Parametrisation
tauMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)*cMax;
tauStep = min(grid.dx, grid.dy)*cMax;

tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;
tStep = dt;

% Sources locations
xS = cat(3, grid.Nx*grid.dx/4, grid.Ny*grid.dy/4);

% Amplitude comparison
grid.newSource(xS, pi/2, pi/2+0.03, nRays, tStep, tMax);
grid.forward_trajectories(1);
grid.forward_amplitudeODE(1);

grid.newSource(xS, pi/2, pi/2+0.03, nRays, tStep, tMax);
grid.forward_trajectories(2);
grid.forward_amplitude(2);
%%  % Sources
%%  for n = 1:nSources
%%      grid.newSource(x{n}, -pi/2, pi/2, nRays, tStep, tMax);
%%      grid.computeHamil(n);
%%      %grid.deleteRays(n);
%%  end

%%  grid.computeForwardParallel(xL, -pi/2, pi/2, nRays, tStep, tMax, true);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex37_gaussBeams;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex37_gaussBeams;

%load gridRT.mat;
%load sensor_data.mat;

%==============================
% Ray trajectories
%==============================
grid.plot_rays(1, 1);

%==============================
% Gaussian beam matrices
%==============================
%grid.plot_YNM(1, 1);
h1 = grid.plot_amplitude(1, 100, true);
figure(h1);
hold on;
semilogy(grid.source(2).phi(1, :), real(grid.source(2).amplitude(1, :)), 'Color', 'r', 'LineWidth', 2);
saveas(gcf, 'Example37_amplitude', 'png');
saveas(gcf, 'Example37_amplitude.fig');

%==============================
% Plot amplitudes
%==============================
nRays = 1;
figure;
semilogy(grid.source(1).phi(1, :), real(grid.source(1).amplitude(1, :)), 'Color', 'b');
hold on;
semilogy(grid.source(2).phi(1, :), real(grid.source(2).amplitude(1, :)), 'Color', 'r', 'LineWidth', 2);

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

%=======================
% Get matrices functions
%=======================
%%  getY = @(index)reshape(permute(grid.source(1).Y(:,index,:),[3,2,1]),3,3);
%%  getN = @(index)reshape(permute(grid.source(1).N(:,index,:),[3,2,1]),3,3);
%%  getM = @(index)reshape(permute(grid.source(1).M(:,index,:),[3,2,1]),3,3);
%%  getP = @(index) cat(1, 1, 1350*permute(grid.source(1).p(:, index, :), [3 2 1]));
%%  defY = @(s, a, b, c) eye(3) + s*[-1 0 0; 0 c^2 0; 0 0 c^2]*[b*i 0 -b*i; 0 a*i 0; -b*i 0 b*i];
%%  defM = @(s, a, b, c) getM(1)*inv(defY(s, a, b, c));

getP = @(index, ray) permute(grid.source(1).p(ray, index, :), [3, 2, 1]);
getDx = @(index, ray) reshape(permute(grid.source(1).Dx(ray,index,:),[3,2,1]), 2, 2);
getDp = @(index, ray) reshape(permute(grid.source(1).Dp(ray,index,:),[3,2,1]), 2, 2);
getD2phase = @(index, ray) reshape(permute(grid.source(1).D2phase(ray,index,:),[3,2,1]), 2, 2);

