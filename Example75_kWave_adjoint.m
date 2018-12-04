% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex75_3D_thinveins_het;

clear all;
close all;

run colourMap;
load ./input_data/sensor_data_veins_14400sensors.mat;

%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
%==============================
% Adjoint
%==============================
% Read results
sensor_data = h5read('output_data/Example75_forward_output_14400sensors.h5', '/p');
save input_data/sensor_data_14400sensors.mat sensor_data;
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
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example75_adjoint_input_14400sensors.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example75_adjoint_input_14400sensors.h5 -o output_data/Example75_adjoint_output_14400sensors.h5 --p_final');


%=========================================================================
% SAVE
%=========================================================================
PML_SIZE = 10;
u_PML = h5read('output_data/Example75_adjoint_output_14400sensors.h5', '/p_final');
u_k_neg = u_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
u_k = max(0, u_k_neg);
% Save pressure
u_matrix = cube2matrix(u_k_neg);
dlmwrite('input_data/pressure_adjoint_kWave_14400sensors.dat', u_matrix, 'delimiter', ' ');

