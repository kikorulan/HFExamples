%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for gridRT class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Measure computational time
tic;
start_time = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 1e3;
dx = 1e-3;
Ny = 1.5e3;
dy = 1e-3;
kernelSize = 100;

grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
v1 = 0.1;
v2 = -0.1;
% Auxiliary matrices
M1 = c0*ones(Nx, floor(Ny/3));
M2 = v1*c0*ones(floor(Nx/2), floor(Ny/3));
M3 = v2*c0*ones(floor(Nx/2), floor(Ny/3));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [M1...
     [M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 

% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;

grid.setCMatrix(cConv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays
nSources = 11;


tauMax = 2*max(Nx*dx, Ny*dy)*c0;
step = 2*min(dx, dy)*c0;

% Alternative shoot rays
% x0 = [[Nx*dx/2]; [Ny*dy*0]];
%% Shoot rays
%for i = 1:1:nRays
%    theta = pi*(i-1)/nRays;
%    p0 = [cos(theta); sin(theta)];
%    grid.computeRay(x0, p0, step, tauMax);
%end

p0 = [0; 1/c0];
%% Shoot rays
for i = 1:1:nSources
    x0 = [Nx*dx*i/nSources-dx; Ny*dy/10];
    grid.newSource(x0, step, tauMax);
    grid.computeRay(i, p0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
for i = 1:1:nSources
    x = grid.getX(i, 1);
    plot(x(1, :), x(2, :));
end
saveas(gcf, 'HighFreq/Examples/Ex08_RT/Example8_Ray_v2.fig');

figure;
x = 1:1:Nx;
y = 1:1:Ny;
flipGrid = grid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, 'HighFreq/Examples/Ex08_RT/Example8_C_v2.fig');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
