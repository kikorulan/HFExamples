%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex62_recon2D_subsample;

close all;
%clear all;

load sensor_data.mat;

run colourMap;
%========================================
% Rgrid definition
%========================================
Nx = 128;           % number of Rgrid points in the x (row) direction
Ny = 256;           % number of Rgrid points in the y (column) direction
dx = 2e-4;        % Rgrid point spacing in the x direction [m]
dy = 2e-4;        % Rgrid point spacing in the y direction [m]
Rgrid = gridRT(Nx, dx, Ny, dy);

% Build domain
Rgrid.setCMatrix(medium.sound_speed);

% Set initial pressure
Rgrid.setUMatrix(sourceKW.p0);

%========================================
% Impulse Response
%========================================
%%  % Set time
%%  dt = 5e-8;
%%  tMax = 4e-5;
%%  Rgrid.setTime(dt, tMax);
%%  % Compute impulse response
%%  Rgrid.impulse_additive('IV');
%%  
%%  save gridRT_impulse.mat Rgrid;

%========================================
% Ray Shooting Parameters
%========================================
load gridRT_impulse;
cMax = max(Rgrid.c(:));
dt = 5e-8;
%dt = min(Rgrid.dx, Rgrid.dy)/cMax/2;

% Measure computational time
tic;
start_time = clock;

% Number of rays & sources
nRays = 2000;% 800

% Parametrisation
tMax = sqrt((Rgrid.Nx*Rgrid.dx)^2+(Rgrid.Ny*Rgrid.dy)^2)/cMax;
tStep = dt;

%========================================
% Forward Problem
%========================================
% Sources locations
subsample_factor = 16;
lenX = 0;
% y = 1
clear x0;
for i = 1:subsample_factor:Nx
    lenX = lenX + 1;
    x0{lenX} = cat(3, (i-1)*Rgrid.dx, 0);
end
% x = 1 & x = end
x_1 = (subsample_factor+1):subsample_factor:Ny;
x_end = fliplr((Ny-subsample_factor):-subsample_factor:1);
for i = 1:length(x_1);
    lenX = lenX + 1;
    if (x_1(i) <= x_end(i))
        x0{lenX}   = cat(3, 0, ((x_1(i)-1)*Rgrid.dy));
        x0{lenX+1} = cat(3, (Rgrid.Nx-1)*dx, ((x_end(i)-1)*Rgrid.dy));
    else
        x0{lenX}   = cat(3, (Rgrid.Nx-1)*dx, ((x_end(i)-1)*Rgrid.dy));
        x0{lenX+1} = cat(3, 0, ((x_1(i)-1)*Rgrid.dy));
    end
    lenX = lenX + 1;
end
% y = end
for i = fliplr(Nx:-subsample_factor:1)
    lenX = lenX + 1;
    x0{lenX} = cat(3, (i-1)*Rgrid.dx, (Rgrid.Ny-1)*Rgrid.dy);
end
nSources = length(x0);


source = Rgrid.computeForwardParallel(x0, 0, 2*pi-0.01, nRays, tStep, tMax, false);
aForward_RT = zeros(nSources, length(source(1).aForward));
for n = 1:nSources
    aForward_RT(n, :) = source(n).aForward;
end

% Save results
%save gridRT.mat Rgrid nRays source x -v7.3;

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex62_recon2D_subsample;

%load gridRT.mat;
%load sensor_data.mat;


%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save RgridRT.mat Rgrid nRays nSources x -v7.3;



