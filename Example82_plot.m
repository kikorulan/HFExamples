% Read data from files
%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex82_TVanalysis;
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex82_TVanalysis;

clear all;
close all;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('./input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%========================================================================================================================
% INITIAL PRESSURE
%========================================================================================================================
u0Matrix = importdata('./input_data/initial_pressure_veins_80x240x240.dat', ' ', 0);
u0 = matrix2cube(u0Matrix, Nz);
h = plot_projection(u0, dx);
a = axes;
t = title('Initial Pressure');
a.Visible = 'off'; 
t.Visible = 'on'; 



%========================================================================================================================
% TV DENOISING
%========================================================================================================================
lambda = 1e-2;
para.maxIter = 200;


% Iter 1
u1 = conTVdenoising(u0, lambda, para);
plot_projection(u1, dx);
colorbar();

% Iter 2
u2 = conTVdenoising(u1, lambda, para);
plot_projection(u2, dx);
colorbar();

% Iter 3
u3 = conTVdenoising(u2, lambda, para);
plot_projection(u3, dx);
colorbar();


% Iter 1 - 600 iters
para_600.maxIter = 600;
u1_600 = conTVdenoising(u0, lambda, para_600);
plot_projection(u1_600, dx);


%========================================================================================================================
% CPP TV DENOISING
%========================================================================================================================
uMatrix = importdata('./output_data/u1.dat', ' ', 0);
u1_cpp = matrix2cube(uMatrix, Nz);
h = plot_projection(u1_cpp, dx);

uMatrix = importdata('./output_data/u2.dat', ' ', 0);
u2_cpp = matrix2cube(uMatrix, Nz);
h = plot_projection(u2_cpp, dx);

uMatrix = importdata('./output_data/u3.dat', ' ', 0);
u3_cpp = matrix2cube(uMatrix, Nz);
h = plot_projection(u3_cpp, dx);


%========================================================================================================================
% TV Differences
%========================================================================================================================
% Maxim
disp('Maximum of u0: ')
maxU0 = max(abs(u0(:)))
disp('Maximum of u1: ')
maxU1 = max(abs(u1(:)))
disp('Maximum of u1_cpp: ')
maxU1_cpp = max(abs(u1_cpp(:)))
disp('Maximum difference u0 and u1')
difU01 = max(abs(u0(:) - u1(:)))
disp('Maximum difference u0 and u1_cpp')
difU01 = max(abs(u0(:) - u1_cpp(:)))
disp('Maximum difference u1 and u1_cpp')
difU11 = max(abs(u1(:) - u1_cpp(:)))

disp('===============================')
% Maxim
disp('Maximum of u1: ')
maxU1 = max(abs(u1(:)))
disp('Maximum of u2: ')
maxU2 = max(abs(u2(:)))
disp('Maximum of u2_cpp: ')
maxU2_cpp = max(abs(u2_cpp(:)))
disp('Maximum difference u1 and u2')
difU12 = max(abs(u1(:) - u2(:)))
disp('Maximum difference u1 and u2_cpp')
difU12_cpp = max(abs(u1(:) - u2_cpp(:)))
disp('Maximum difference u2 and u2_cpp')
difU22 = max(abs(u2(:) - u2_cpp(:)))

disp('===============================')
% Maxim
disp('Maximum of u2: ')
maxU2 = max(abs(u2(:)))
disp('Maximum of u3: ')
maxU3 = max(abs(u3(:)))
disp('Maximum of u3_cpp: ')
maxU3_cpp = max(abs(u3_cpp(:)))
disp('Maximum difference u2 and u3')
difU23 = max(abs(u2(:) - u3(:)))
disp('Maximum difference u2 and u3_cpp')
difU23_cpp = max(abs(u2(:) - u3_cpp(:)))
disp('Maximum difference u3 and u3_cpp')
difU33 = max(abs(u3(:) - u3_cpp(:)))


%========================================================================================================================
% TV convergence
%========================================================================================================================
lambda = '5e-4';
iter = 99;

uMatrix = importdata(['./results/u0.dat'], ' ', 0);
u_ini = matrix2cube(uMatrix, Nz);
h = plot_projection(u_ini, dx);

uMatrix = importdata(['./results/L', lambda, '/u', int2str(iter), '_lambda', lambda, '.dat'], ' ', 0);
u_TV = matrix2cube(uMatrix, Nz);
h = plot_projection(u_TV, dx);

uMatrix = importdata(['./results/L', lambda, '/u', int2str(iter+1), '_lambda', lambda, '.dat'], ' ', 0);
u_TVnext = matrix2cube(uMatrix, Nz);
h = plot_projection(u_TVnext, dx);


dfIni_iter = max(abs(u_ini(:) - u_TV(:)))
dfIni_last = max(abs(u_ini(:) - u_TVnext(:)))
dfIter = max(abs(u_TV(:) - u_TVnext(:)))

