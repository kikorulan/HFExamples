%cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex81_3D_veins_subsampled_het;
cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex81_3D_veins_subsampled_het;

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


%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('./output_data/Example81_forward_output_3600sensors.h5', '/p');

%==================================================
% Plot
%==================================================
x = importdata('./output_data/x1_subiter0_iter1.dat', ' ', 0);
x1 = matrix2cube(x, Nz);
x = importdata('./output_data/x2_subiter0_iter1.dat', ' ', 0);
x2 = matrix2cube(x, Nz);
x = importdata('./output_data/x3_subiter0_iter1.dat', ' ', 0);
x3 = matrix2cube(x, Nz);

plot_projection(x1, 1);
plot_projection(x2, 1);
plot_projection(x3, 1);


y0 = importdata('./output_data/y0_subiter0_iter1.dat', ' ', 0);
y1 = importdata('./output_data/y1_subiter0_iter1.dat', ' ', 0);
y2 = importdata('./output_data/y2_subiter0_iter1.dat', ' ', 0);
y3 = importdata('./output_data/y3_subiter0_iter1.dat', ' ', 0);
y4 = importdata('./output_data/y4_subiter0_iter1.dat', ' ', 0);
y5 = importdata('./output_data/y5_subiter0_iter1.dat', ' ', 0);


figure;
plot(y0);
hold on;
plot(y1);
plot(y2);
plot(y3);
plot(y4);
plot(y5);
legend('y0', 'y1', 'y2', 'y3', 'y4', 'y5')


index = 1;
y0_kWave = sensor_data(index, :);
figure;
plot(y0/max(y0));
hold on;
plot(y0_kWave/max(y0_kWave));
% Forward signal
time_signal = importdata(['./input_data/forwardSignal_reference_3600sensors.dat'], ' ', 0);
y0_ref = time_signal(2:end, :);
figure
plot(y0_ref(1, :))
