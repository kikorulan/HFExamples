% Example for gridFFM class

clear all
close all

% Measure computational time
tic;
time1 = clock;

%  % Grid definition
Nx = 200;
dx = 1e-3;
Ny = 200;
dy = 1e-3;

grid = gridFMM(Nx, dx, Ny, dy);
grid.setCMatrix(1500*ones(Nx, Nx));

% Initial conditions and run grid
v = [floor(Nx/2) floor(2*Ny/3) 1e-6; floor(Nx/3) floor(Ny/3) 1e-6; floor(2*Nx/3) floor(Ny/3) 1e-6];
grid.setInitialConditions(v);
grid.runFMM();

% Measure computational time
time2 = clock;
disp(['  Grid computation time ' num2str(etime(time2, time1))]);

% Compute Error
gradPhase = NaN(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        gradPhase(i, j) = norm(grid.gradPhase([i j]));
    end
end
gradPhase0 = 1./grid.c;
errorGradPhase = abs(gradPhase0 - gradPhase)./gradPhase0;


time3 = clock;
disp(['  Error computation time ' num2str(etime(time3, time2))]);

%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%
x = 1:1:Nx;
y = 1:1:Ny;

% Phase
%  grid.setU([v(1, 1) v(1, 2)], NaN);
%  grid.setU([v(2, 1) v(2, 2)], NaN);
%  grid.setU([v(3, 1) v(3, 2)], NaN);

figure;
surf(x, y, grid.u, 'EdgeColor', 'none');
view(2);
saveas(gcf, 'Example2_Phase.fig');

% Error
%  errorGradPhase(v(1, 1), v(1, 2)) = NaN;
%  errorGradPhase(v(2, 1), v(2, 2)) = NaN;
%  errorGradPhase(v(3, 1), v(3, 2)) = NaN;

figure;
surf(x, y, errorGradPhase, 'EdgeColor', 'none');
view(2);
saveas(gcf, 'Example2_Error.fig');
maxErrorGradPhase = max(max(errorGradPhase))
