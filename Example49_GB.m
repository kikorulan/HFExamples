%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex49_GB;

close all;
clear all;

run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 257;           % number of grid points in the y (column) direction
dx = 1e-3;        % grid point spacing in the x direction [m]
dy = 1e-3;        % grid point spacing in the y direction [m]

% Sound speed
c0 = 1500;
w0 = 1e8;

% Initialize gaussian beams
sVector = (0:dx:(Nx-1)*dx)';
sMatrix = repmat(sVector, [1 Ny]);

% Q & P
p0 = 1;
q0 = i*w0*dx*dx/2;
qMatrix = c0*p0*sMatrix + q0;

% Distance matrix
yVector = (1:1:Ny);
yHalf = floor(Ny/2) + 1;
dVector = abs(yVector - yHalf)*dx;
dMatrix = repmat(dVector, [Nx 1]);

% GB matrix
GBmatrix = exp(-i*w0*(0*sMatrix + p0/2./qMatrix.*dMatrix.*dMatrix));

%================================================================================
% Plot results
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex49_GB;
position = [700 700 400 600];

% S matrix
figure;
imagesc(sMatrix);
colorbar();
title('S')

% D matrix
figure;
imagesc(dMatrix);
colorbar();
title('D')

% GB matrix
GBmatrixR = abs(GBmatrix);
figure;
imagesc(GBmatrixR);
colorbar();
%dim = [.1 .1 .3 .3];
%str = 'Straight Line Plot from 1 to 10';
%annotation('textbox',dim,'String',str,'FitBoxToText','on', 'Color', [1 1 1]);
set(gcf, 'pos', position);
title({'GB', 'p0 = 1, q0 = i*w0*dx*dx/2', 'w0 = 1e8'})
saveas(gcf, 'Example49_1e8', 'epsc');

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
