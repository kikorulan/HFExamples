% Example for gridFFM class

clear all
close all

% Measure computational time
tic;
t1 = clock;

% Grid definition
Nx = 300;
dx = 1e-3;
Ny = 300;
dy = 1e-3;
c = 1500;
u0 = 1e-5;
grid = gridFMM(Nx, dx, Ny, dy);
grid.setCMatrix(c*ones(Nx, Ny));

% Initial conditions and run grid
v = [floor(Nx/2) floor(2*Ny/3) [u0-0.8e-5]; floor(Nx/3) floor(Ny/3) [u0-0.4e-5]; floor(2*Nx/3) floor(Ny/3) u0];
grid.setInitialConditions(v);
grid.runFMM();

%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%
% Axis definition
x = 1:1:Nx;
y = 1:1:Ny;

h = figure;
surf(x, y, grid.u, 'EdgeColor', 'none');
view(2);
saveas(h, 'Example4_Forward.fig');

t2 = clock;
% Measure computational time for forward problem
disp(['  total computation time ' num2str(etime(t2, t1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ISOSURFACE AS INITIAL CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2 = clock;
% Find isosurface
incTime = 1.1*sqrt(dx*dx + dy*dy)/c;
T0 = 6e-5;
M = grid.findIsosurface(T0, incTime);
h = figure;
x = 1:1:Nx;
y = 1:1:Ny;
surf(x, y, M, 'EdgeColor', 'none');
view(2);
saveas(h, 'Example4_IsoSurf.fig');

% Consider the isosurface as the initial condition
gridIso = gridFMM(Nx, dx, Ny, dy);
gridIso.setCMatrix(c*ones(Nx, Ny));
gridIso.setInitialConditionsMatrix(M);


% Run in 3 steps: each step of 1e-5 s
    % First step: 1e-5
    gridIso.runFMM_Time(1e-5);
    h = figure;
    surf(x, y, gridIso.u, 'EdgeColor', 'none');
    view(2);
    colorbar;
    saveas(gcf, 'Example4_Backward-1.fig');

    % Second step: 2e-5
    gridIso.runFMM_Time(3e-5);
    h = figure;
    surf(x, y, gridIso.u, 'EdgeColor', 'none');
    view(2);
    colorbar;
    saveas(gcf, 'Example4_Backward-2.fig');

    % Third step: 3e-5
    gridIso.runFMM_Time(6e-5);
    h = figure;
    surf(x, y, gridIso.u, 'EdgeColor', 'none');
    view(2);
    colorbar;
    saveas(gcf, 'Example4_Backward-3.fig');

% Plot error
gridIso_inv = -gridIso.u + max(gridIso.u(gridIso.u~=Inf)) + u0;
gridError = log10(abs(grid.u - gridIso_inv)/max(grid.u(grid.u~=inf)));
h = figure;
surf(x, y, gridError, 'EdgeColor', 'none');
view(2);
colorbar;
saveas(gcf, 'Example4_Error.fig');
% 

t3 = clock;
% Measure computational time for backward problem
disp(['  total computation time ' num2str(etime(t3, t2))]);
