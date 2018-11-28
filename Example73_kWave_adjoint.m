% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex73_3D_heterogeneous_sync;

clear all;
close all;

run colourMap;
load ./input_data/sensor_data_veins_25sensors.mat;

%=========================================================================
% SIMULATION
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';
%==============================
% Adjoint
%==============================
nS = 13;
% Read results
sensor_data = h5read('output_data/Example73_forward_output_25sensors.h5', '/p');
save input_data/sensor_data_25sensors.mat sensor_data;
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);
index = sensor_index(nS);
% Consider all sensors
source_adjoint.p_mask = sensor.mask;
source_adjoint.p_mask(:) = 0;
source_adjoint.p_mask(index) = 1;
source_adjoint.p = fliplr(sensor_data(nS, :));
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example73_adjoint_input_25sensors.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example73_adjoint_input_25sensors.h5 -o output_data/Example73_adjoint_output_25sensors.h5 --p_final');

