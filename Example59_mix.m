% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex59_3D_4balls;

clear all;
close all;

run colourMap;
load input_data/sensor_data_4balls.mat;

%=========================================================================
% k-Wave
%=========================================================================
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/GCC-5.4/lib64';

% Read results - RT data
sensor_data_RT = importdata('output_data/ForwardSignal_2168sensors.dat', ' ', 0);
% Number of sensors
sensor_index = find(sensor.mask == 1);
nSensors = length(sensor_index);

% Consider all sensors
source_adjoint.p_mask = sensor.mask;
source_adjoint.p = fliplr(sensor_data_RT(2:end, :));
% Sensor
sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor_adjoint.record = {'p_final'};
% Save and run
kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example59_adjoint_input.h5');
system('../kspaceFirstOrder3D-OMP -i input_data/Example59_adjoint_input.h5 -o output_data/Example59_adjoint_output_RTdata.h5 --p_final');


%=========================================================================
% Ray Tracing
%=========================================================================
% Write results
sensor_data_kWave = h5read('output_data/Example59_forward_output.h5', '/p');
sensor_data_kWave = [sensor_data_RT(1, :); sensor_data_kWave];
dlmwrite('input_data/forwardSignal_2168sensors_kWave.dat', sensor_data_kWave, 'delimiter', ' ');



