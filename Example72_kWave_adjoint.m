% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex72_3D_veins_heterogeneous;

clear all;
close all;

run colourMap;
load ./input_data/sensor_data_veins_1600sensors.mat;

%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
%==============================
% Adjoint
%==============================
% Read results
sensor_data = h5read('output_data/Example72_forward_output_1600sensors.h5', '/p');
save input_data/sensor_data_1600sensors.mat sensor_data;
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
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example72_adjoint_input_1600sensors.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example72_adjoint_input_1600sensors.h5 -o output_data/Example72_adjoint_output_1600sensors.h5 --p_final');

