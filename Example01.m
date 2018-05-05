% Example for gridFFM class

clear all
close all

% Measure computational time
tic;
start_time = clock;

% Grid definition
Nx = 50;
dx = 0.1;
Ny = 50;
dy = 0.1;
grid = gridFMM(Nx, dx, Ny, dy);
v = [25 25 1];
grid.setInitialConditions(v);

grid.runFMM();

% Plot results
figure;
Z = grid.phase;
x = 1:1:Nx;
y = 1:1:Ny;
surf(x, y, Z);

view(2);
saveas(gcf, 'Example1', 'epsc');

end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
