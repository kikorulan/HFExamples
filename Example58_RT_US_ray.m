% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex58_US;

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
factor = 0.1;
% Build Peaks
frame = 2.5;
x = -frame:(2*frame)/(Ny-1):frame;
y = -frame:(2*frame)/(Nx-1):frame;
[X, Y] = meshgrid(x, y);
p = peaks(X, Y);
%cMatrix = c0*(ones(Nx, Ny) + factor*p/max(p(:)));
cMatrix = c0*ones(Nx, Ny);
grid.setCMatrix(cMatrix);

%========================================
% Impulse Response
%========================================
% Set time
cMax = max(grid.c(:));
cMin = min(grid.c(:));
dt = min(grid.dx, grid.dy)/cMax/3;
tMax = 3e-5;
grid.setTime(dt, tMax);

%=======================================
% Create Sources
%=======================================
% Number of rays & sources
nRays = 2000; % 3000
nSources = 1;
% Maximum t and step size
tMax = sqrt((grid.Nx*grid.dx)^2+(grid.Ny*grid.dy)^2)/cMin;
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

%=======================================
% Compute Lengths
%=======================================
% Create the new sources
prop_time_vec = [];
for i = 1:nSources
    disp(strcat('Estimation ', int2str(i)));
    [index_vec, prop_time] = grid.ssestimation_rays(source, i);
    prop_time_vec = [prop_time_vec; prop_time];
end
grid.ssestimation_allocate(source);
for i = 1:nSources
    disp(i);
    grid.ssestimation_length(source, i);
end


%=======================================
% Verify computation
%=======================================
%%  %load prop_time_126;
%%  load prop_time_homogeneous_126;
%%  
%%  %prop_time_vec = reshape(prop_time_matrix, [], 1);
%%  prop_time_vec_homo = reshape(prop_time_matrix_homo, [], 1);
%%  c_vec = reshape(grid.c, [], 1);
%%  n_vec = 1./c_vec;
%%  
%%  % Difference with respect to heterogeneous
%%  dif_hete = abs(grid.LMatrix*n_vec - prop_time_vec);
%%  maxDifHete = max(dif_hete)
%%  normDifHete = norm(dif_hete)
%%  % Difference with respect to homogeneous
%%  dif_homo = abs(grid.LMatrix*n_vec - prop_time_vec_homo);
%%  maxDifHomo = max(dif_homo)
%%  normDifHomo = norm(dif_homo)

%=========================================================================================
% VISUALIZATION
%=========================================================================================
%=======================================
% Plot rays
%=======================================
source(1).plot_rays(grid, 300);

%=======================================
% Comparison between time estimation
%=======================================
%%  % Propagation times
%%  figure;
%%  plot(prop_time_vec, 'Color', 'r');
%%  hold on;
%%  plot(prop_time_vec_homo, 'Color', 'b');
%%  
%%  % Difference in estimation
%%  figure;
%%  plot(dif_hete, 'Color', 'r');
%%  hold on;
%%  plot(dif_homo, 'Color', 'b');
%%  legend('Heterogeneous', 'Homogeneous');

%=======================
% Time
%=======================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

