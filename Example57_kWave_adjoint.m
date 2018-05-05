% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex57_RT_3D_5x5x5;

clear all;
close all;

run colourMap;
load input_data/sensor_data_3balls.mat;

%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/GCC-5.4/lib64';
%==============================
% Adjoint
%==============================
% Read results
sensor_data = h5read('output_data/Example57_forward_output.h5', '/p');
save input_data/sensor_data.mat sensor_data;
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);

% Consider all sensors
source_adjoint.p_mask = sensor.mask;
source_adjoint.p = fliplr(sensor_data);
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example57_adjoint_input.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example57_adjoint_input.h5 -o output_data/Example57_adjoint_output.h5 --p_final');

%=========================================================================
% VISUALISATION
%=========================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex57_RT_3D_5x5x5;

% Axis
x_axis = 0:kgrid.dx:(kgrid.Nx-1)*kgrid.dx;
y_axis = 0:kgrid.dy:(kgrid.Ny-1)*kgrid.dy;

%==============================
% Reconstruction
%==============================
p0_recon_PML = h5read('output_data/Example57_adjoint_output.h5', '/p_final');
PML_size = 10;
p0_recon = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));

% Normalisation - kWave
inputKWave = p0_recon(:, :, 64);
normKWave = max(inputKWave(:));
inputKWave_norm = inputKWave/normKWave;

% Plot figure
figure;
surf(x_axis, y_axis, inputKWave_norm, 'EdgeColor', 'none');
axis([0 x_axis(end) 0 y_axis(end)]);
view(2);
colorbar();
xlabel('y (m)');
ylabel('x (m)');
title('kWave recon - z = 0.063');
saveas(gcf, 'Example57_kWave_recon', 'png');
saveas(gcf, 'Example57_kWave_recon.fig');

cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples;
