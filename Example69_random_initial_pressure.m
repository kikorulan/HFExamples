% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex69_3D_RT_norm;

clear all;
close all;

buildDomain = 0;
plotPowerIt = 1;

%=========================================================================
% DOMAIN DEFINITION
%=========================================================================
Nx = 80;           % number of grid points in the x (row) direction
Ny = 240;          % number of grid points in the y (column) direction
Nz = 240;          % number of grid points in the z (column) direction
boundary = 10;
dx = 5.3e-5;

if (buildDomain)
rng(1); % Set random generator
initial_pressure = zeros(Nx, Ny, Nz);
initial_pressure(boundary+1:Nx-boundary, boundary+1:Ny-boundary, boundary+1:Nz-boundary) = rand(Nx - 2*boundary, Ny - 2*boundary, Nz - 2*boundary);

plot_projection(initial_pressure, 1);

% Save initial pressure
initial_pressure_matrix = cube2matrix(initial_pressure);
dlmwrite('input_data/initial_pressure_random.dat', initial_pressure_matrix, 'delimiter', ' ');
% Sound speed
c0 = 1580.0001;
c = c0*ones(Nx, Ny, Nz);
c_matrix = cube2matrix(c);
dlmwrite('input_data/sound_speed.dat', c_matrix, 'delimiter', ' ');
end


%=========================================================================
% PLOT POWER IT
%=========================================================================
if (plotPowerIt)

dot_product = @(x, y) sqrt(sum(x(:).*y(:)));

% Load initial pressure
initial_pressure_mat = importdata('input_data/initial_pressure_random.dat', ' ', 0);
initial_pressure = max(0, matrix2cube(initial_pressure_mat, Nz));
dot_product_ip = dot_product(initial_pressure, initial_pressure);

% Adjoint
pressure_adjoint_mat = importdata('input_data/pixelPressure_initial.dat', ' ', 0);
pressure_adjoint = max(0, matrix2cube(pressure_adjoint_mat, Nz));
norm_initial_pow = dot_product(initial_pressure, pressure_adjoint)/dot_product_ip;
norm_initial = nthroot(norm_initial_pow, 1)

% ITER 0
pressure_adjoint_mat = importdata('input_data/pixelPressure_iter0.dat', ' ', 0);
pressure_adjoint = max(0, matrix2cube(pressure_adjoint_mat, Nz));
norm_pow = dot_product(initial_pressure, pressure_adjoint)/dot_product_ip;
norm_iter0 = nthroot(norm_pow, 2)

% ITER 1
pressure_adjoint_mat = importdata('input_data/pixelPressure_iter1.dat', ' ', 0);
pressure_adjoint = max(0, matrix2cube(pressure_adjoint_mat, Nz));
norm_pow = dot_product(initial_pressure, pressure_adjoint)/dot_product_ip;
norm_iter1 = nthroot(norm_pow, 3)

% ITER 2
pressure_adjoint_mat = importdata('input_data/pixelPressure_iter2.dat', ' ', 0);
pressure_adjoint = max(0, matrix2cube(pressure_adjoint_mat, Nz));
norm_pow = dot_product(initial_pressure, pressure_adjoint)/dot_product_ip;
norm_iter2 = nthroot(norm_pow, 3)

% ITER 3
pressure_adjoint_mat = importdata('input_data/pixelPressure_iter3.dat', ' ', 0);
pressure_adjoint = max(0, matrix2cube(pressure_adjoint_mat, Nz));
norm_pow = dot_product(initial_pressure, pressure_adjoint)/dot_product_ip;
norm_iter3 = nthroot(norm_pow, 4)

% ITER 4
pressure_adjoint_mat = importdata('input_data/pixelPressure_iter4.dat', ' ', 0);
pressure_adjoint = max(0, matrix2cube(pressure_adjoint_mat, Nz));
norm_pow = dot_product(initial_pressure, pressure_adjoint)/dot_product_ip;
norm_iter4 = nthroot(norm_pow, 5)



end
