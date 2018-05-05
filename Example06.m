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
Nx = 100;
dx = 1e-3;
Ny = 100;
dy = 1e-3;
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
% Auxiliary matrices
c = c0*ones(Nx, Ny);

grid.setCMatrix(c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays
nRays = 30;


tauMax = Nx*dx+Ny*dy;
step = 5e-4;

%%%%%%
x0 = [[Nx*dx/3]; [2*Ny*dy/3]]; % Ray 1
x1 = [[2*Nx*dx/3]; [Ny*dy/3]]; % Ray 2
% Shoot rays
for i = 1:1:nRays
    theta = 2*pi*(i-1)/nRays;
    p0 = [cos(theta); sin(theta)];
    grid.computeRay(x0, p0, step, tauMax);
    grid.computeRay(x1, p0, step, tauMax);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
for i = 1:2:2*nRays
    plot(grid.ray(i).x(1, :), grid.ray(i).x(2, :), 'r');
    plot(grid.ray(i+1).x(1, :), grid.ray(i+1).x(2, :));
end
saveas(gcf, 'HighFreq/Examples/Ex6/Example6_Ray.fig');

figure;
x = 1:1:Nx;
y = 1:1:Ny;
flipGrid = grid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, 'HighFreq/Examples/Ex6/Example6_C.fig');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
