% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex45_photo;

clear all;
close all;

run colourMap;
%=========================================================================
% DOMAIN DEFINITION
%=========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
Nz = 128;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]
dz = 1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 2000;
v1 = -.5;
radi = Nx/5;
kernelSize = 30;
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
dimZ = Nz + 2*kernelSize;
c = c0*ones(dimX, dimY, dimZ);
c = addSphere(c,   floor(dimX/2),   floor(dimY/2), floor(dimZ/2), radi, c0*v1);

cSlice = c(:, :, floor(dimZ/2));
K = ones(kernelSize, kernelSize, kernelSize);
cConvSlice = convn(cSlice, K, 'same')/kernelSize/kernelSize/kernelSize;
cConvSlice = cConvSlice(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
cConv = convn(c, K, 'same')/kernelSize/kernelSize/kernelSize;
cConv = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);

% Build initial pressure
%c_matrix = cube2matrix(repmat(cConvSlice, [1 1 Nz]));
c_matrix = cube2matrix(cConv);
dlmwrite('soundSpeed_ball.dat', c_matrix, 'delimiter', ' ');

%=========================================================================
% SIMULATION
%=========================================================================
%% Call C++ code
%setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/GCC-5.4/lib64';

system('./RTsolver dimensions.dat soundSpeed_ball.dat initialPressure_3balls.dat sensors.dat');

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex45_photo;


%==================================================
% TRAJECTORIES
%==================================================
% Import data
filenameData = 'output_data/Trajectory0.dat';
trajectories = importdata(filenameData, ' ', 0);

% Read number of rays and steps
[nSteps nRays] = size(trajectories);
%nRays = floor(nRays/2);
xCoord = trajectories(:, 1:3:nRays);
yCoord = trajectories(:, 2:3:nRays);
zCoord = trajectories(:, 3:3:nRays);

nRays = floor(nRays/3);
% Plot the figure
figure;
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Nx*dx 0 Ny*dy 0 Nz*dz]);
hold on;
colours = winter(nRays);
for n = 1:nRays
    plot3(xCoord(:, n), yCoord(:, n), zCoord(:, n), 'Color', colours(n, :));
end

% Select curved trajectories
curveTrajecX = xCoord(2:end, :) - xCoord(1:end-1, :);
curveTrajecY = yCoord(2:end, :) - yCoord(1:end-1, :);
curveTrajecZ = zCoord(2:end, :) - zCoord(1:end-1, :);

binaryX = (abs(max(curveTrajecX)-min(curveTrajecX))./abs(max(curveTrajecX))<0.3);
binaryY = (abs(max(curveTrajecY)-min(curveTrajecY))./abs(max(curveTrajecY))<0.3);
binaryZ = (abs(max(curveTrajecZ)-min(curveTrajecZ))./abs(max(curveTrajecZ))<0.3);

binaryTrajec = ~(binaryX&binaryY&binaryZ);
xCoordSub = xCoord(:, binaryTrajec);
yCoordSub = yCoord(:, binaryTrajec);
zCoordSub = zCoord(:, binaryTrajec);

figure;
ax = gca;
ax.GridAlpha = 1;
grid on;
axis([0 Nx*dx 0 Ny*dy 0 Nz*dz]);
hold on;
nRaysSub = sum(binaryTrajec);
colours = winter(nRaysSub);
for n = 1:nRaysSub
    plot3(xCoordSub(:, n), yCoordSub(:, n), zCoordSub(:, n), 'Color', colours(n, :));
end

%==================================================
% GRADIENT
%==================================================
% Import data
delimiterIn = ' ';
headerlinesIn = 0;
n_sound_speed = importdata('output_data/n.dat', delimiterIn, headerlinesIn);
gradX = importdata('output_data/gradX.dat', delimiterIn, headerlinesIn);
gradY = importdata('output_data/gradY.dat', delimiterIn, headerlinesIn);
gradZ = importdata('output_data/gradZ.dat', delimiterIn, headerlinesIn);

% Plot gradient
figure;
surf(n_sound_speed, 'EdgeColor', 'none');
view(2); 
title('N ');

figure;
surf(gradX, 'EdgeColor', 'none');
view(2);
title('Grad X');

figure;
surf(gradY, 'EdgeColor', 'none');
view(2);
title('Grad Y');

figure;
surf(gradZ, 'EdgeColor', 'none');
view(2);
title('Grad Z');


cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;
