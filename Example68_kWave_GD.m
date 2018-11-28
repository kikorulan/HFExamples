% Heterogeneous Propagation Medium Example
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex68_3D_veins_resize;

clear all;
close all;

% LIBRARIES
setenv LD_LIBRARY_PATH '/cs/research/medim/projects2/projects/frullan/lib/root/lib64';

%====================
% LOAD DATA
%====================
load ./input_data/sensor_data_veins_1600sensors.mat;
% smooth the initial pressure distribution and restore the magnitude
PML_size = 10;
%u_PML = h5read('output_data/Example68_adjoint_output_1600sensors.h5', '/p_final');
u_PML = h5read('output_data/Example68_adjoint_output_57600sensors.h5', '/p_final');
u_k_neg = u_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
u_k = max(0, u_k_neg);


% Save pressure
u_matrix = cube2matrix(u_k_neg);
dlmwrite('input_data/pressure_adjoint_kWave_57600sensors.dat', u_matrix, 'delimiter', ' ');


%============================================================
% SIMULATION
%============================================================

% PARAMETERS
tau = 10;
lambda = 0.005;
nIter = 1;

for n = 1:nIter
    %====================
    % FORWARD
    %====================
    clear source;
    source.p0 = smooth(kgrid, u_k, true);
    source.p0 = max(0, source.p0);
    % Save to disk
    filename = 'input_data/Example68_forward_input_A.h5';
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', filename);
    % Call C++ code
    system('../kspaceFirstOrder3D-OMP -i input_data/Example68_forward_input_A.h5 -o output_data/Example68_forward_output_A.h5');
    
    %==============================
    % ADJOINT
    %==============================
    load ./input_data/sensor_data_1600sensors;
    % Read results
    sensor_data_adjoint = h5read('output_data/Example68_forward_output_A.h5', '/p');
    % Number of sensors
    sensor_index = find(sensor.mask == 1);
    nSensors = length(sensor_index);
    
    % Consider all sensors
    source_adjoint.p_mask = sensor.mask;
    source_adjoint.p = fliplr(sensor_data_adjoint - sensor_data);
    % Sensor
    sensor_adjoint.mask = ones(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    sensor_adjoint.record = {'p_final'};
    % Save and run
    kspaceFirstOrder3D(kgrid, medium, source_adjoint, sensor_adjoint, input_args{:}, 'SaveToDisk', 'input_data/Example68_adjoint_input_A.h5');
    system('../kspaceFirstOrder3D-OMP -i input_data/Example68_adjoint_input_A.h5 -o output_data/Example68_adjoint_output_A.h5 --p_final');
    
    %==============================
    % UPDATE 
    %==============================
    p_PML = h5read('output_data/Example68_adjoint_output_A.h5', '/p_final');
    p = p_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size);
    u_k = softThresh(u_k - 2*tau*p, lambda);
    u_k = max(0, u_k);

    %==============================
    % PLOT
    %==============================
    plot_projection(u_k, kgrid.dx);
    pause(0.1);
end

