% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex59_RT_US_ray;

close all;
%clear all;

% Measure computational time
tic;
start_time = clock;

%=======================================
% Grid definition
%=======================================
Nx = 128;           
Ny = 128;           
dx = 1e-4;        
dy = 1e-4;        
grid = gridRT(Nx, dx, Ny, dy);

% Build domain
c0 = 1500;
cMatrix = 1500*ones(Nx, Ny);
grid.setCMatrix(cMatrix);
save cMatrix.mat cMatrix;

%========================================
% Impulse Response
%========================================
% Set time
cMax = max(grid.c(:));
dt = min(grid.dx, grid.dy)/cMax/3;
tMax = 2e-5;
grid.setTime(dt, tMax);

%=======================================
% Ray Shooting
%=======================================
% Number of rays & sources
nRays = 2000; % 3000
nSources = 1;
% Maximum t and step size
tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMax;
tStep = dt;

% Sources locations
subsampleFactor = 4;
lenX = 0;
% y = 1
clear x0;
for i = 1:subsampleFactor:Nx
    lenX = lenX + 1;
    x0{lenX} = cat(3, (i-1)*grid.dx, 0);
end
% x = 1 & x = end
x_1 = (subsampleFactor+1):subsampleFactor:Ny;
x_end = fliplr((Ny-subsampleFactor):-subsampleFactor:1);
for i = 1:length(x_1);
    lenX = lenX + 1;
    if (x_1(i) <= x_end(i))
        x0{lenX}   = cat(3, 0, ((x_1(i)-1)*grid.dy));
        x0{lenX+1} = cat(3, (grid.Nx-1)*dx, ((x_end(i)-1)*grid.dy));
    else
        x0{lenX}   = cat(3, (grid.Nx-1)*dx, ((x_end(i)-1)*grid.dy));
        x0{lenX+1} = cat(3, 0, ((x_1(i)-1)*grid.dy));
    end
    lenX = lenX + 1;
end
% y = end
for i = fliplr(Nx:-subsampleFactor:1)
    lenX = lenX + 1;
    x0{lenX} = cat(3, (i-1)*grid.dx, (grid.Ny-1)*grid.dy);
end

angleMin = 0;
angleMax = 2*pi-0.01;

nSources = length(x0);
clear source;
for i = 1:nSources
    source(i) = grid.newSource(x0{i}, angleMin, angleMax, nRays, tStep, tMax);
end
% Create the new sources
for i = 1:nSources
    disp(strcat('Estimation ', int2str(i)));
    index_vec = grid.ssestimation_rays(source, i);
end
% Compute the rays and amplitudes
disp(strcat('Source ', int2str(1)));

%=========================================================================================
% VISUALIZATION
%=========================================================================================
%=======================================
% Plot rays
%=======================================
source(1).plot_rays(grid, 300);

%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

