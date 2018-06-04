%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex61_density;

close all;
%clear all;


run colourMap;
%========================================
% Grid definition
%========================================
Nx = 128;           % number of grid points in the x (row) direction
Ny = 256;           % number of grid points in the y (column) direction

%==============================
% define the properties of the propagation medium
%==============================
% Build domain
c0 = 1;
v1 = 10;
v2 = 0.1;
kernelSize = 30;
% Auxiliary matrices
dimX = Nx + 2*kernelSize;
dimY = Ny + 2*kernelSize;
M1 = c0*ones(dimX, dimY);
c = M1;
c = addCircle(c, kernelSize + floor(Nx/2), kernelSize + floor(Ny/2), floor(Nx/6), c0*v1);
% Kernel convolution
K = ones(kernelSize);
cConv = conv2(c, K, 'same')/kernelSize/kernelSize;
medium.sound_speed = cConv(1+kernelSize:end-kernelSize, 1+kernelSize:end-kernelSize);
medium.density = 1;
c = medium.sound_speed;
%========================================
% Ray Shooting
%========================================
%load gridRT_impulse;

% Measure computational time
tic;
start_time = clock;
%integral = density_integral(c, [64; 1]);
integral = density_integral(c);

%==================================================================================
% Plot results
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex46_caustics;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;


%==============================
% Travel times
%==============================
figure;
surf(integral{1}', 'EdgeColor', 'none');
view(2);


%==============================
% Sound speed + rays
%==============================
h = figure;
hold on;
imagesc(c');
xlabel('x (m)');
ylabel('y (m)');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

