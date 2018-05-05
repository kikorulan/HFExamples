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
c0 = 150;
v1 = 0.3;
v2 = -0.3;
% Auxiliary matrices
M1 = c0*ones(Nx, floor(Ny/2));
M2 = v1*c0*ones(floor(Nx/2), floor(Ny/2));
M3 = v2*c0*ones(floor(Nx/2), floor(Ny/2));
M4 = triu(M2);
M5 = fliplr(triu(M3)')';
c = [[M1 + [M4; M5]]...
     [M1 + [M2; M3]]]; 



grid.setCMatrix(c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ray Shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of rays
nRays = 53;


tauMax = Nx*dx+Ny*dy;
step = 5e-4;

% Alternative shoot rays
% x0 = [[Nx*dx/2]; [Ny*dy*0]];
%% Shoot rays
%for i = 1:1:nRays
%    theta = pi*(i-1)/nRays;
%    p0 = [cos(theta); sin(theta)];
%    grid.computeRay(x0, p0, step, tauMax);
%end

p0 = [0; 1];
%% Shoot rays
for i = 1:1:nRays
    x0 = [Nx*dx*i/nRays-dx; 0];
    grid.computeRay(x0, p0, step, tauMax);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
for i = 1:1:nRays
    plot(grid.ray(i).x(1, :), grid.ray(i).x(2, :));
end
saveas(gcf, 'HighFreq/Examples/Ex5/Example5_Ray.fig');

figure;
x = 1:1:Nx;
y = 1:1:Ny;
flipGrid = grid.n';
surf(x, y, flipGrid, 'EdgeColor', 'none');
view(2);
saveas(gcf, 'HighFreq/Examples/Ex5/Example5_C.fig');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
