% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex69_3D_RT_norm;

clear all;
close all;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', ' ', 0);
Nx = dim(1, 1); dx = dim(2, 1);
Ny = dim(1, 2); dy = dim(2, 2);
Nz = dim(1, 3); dz = dim(2, 3);

%============================================================================================================================================
% FORWARD PROBLEM
%============================================================================================================================================
nSensors = 16;

%==================================================
% TIME SIGNAL - kWave
%==================================================
sensor_data = h5read('./output_data/Example69_forward_output_100sensors_homo.h5', '/p');
figure;
imagesc(sensor_data);
%%  hold on;
%%  for n = 1:nSensors
%%      plot(sensor_data(n, :));
%%  end

%==================================================
% TIME SIGNAL - RT
%==================================================
% Import data
filenameData = 'output_data/ForwardSignal.dat';
timeSignal = importdata(filenameData, ' ', 0);
timeRT = timeSignal(1, :);
inputRT = timeSignal(2:end, :);
% Plot
figure;
imagesc(inputRT);
%%  hold on;
%%  for n = 1:nSensors
%%      plot(inputRT(n, :));
%%  end

%==================================================
% Comparison
%==================================================
for n = 1:10
    figure;
    hold on;
    plot(sensor_data(n, :), 'Color', 'r', 'LineWidth', 2);
    plot(inputRT(n, :), 'Color', 'b');
end

%==================================================
% COMPUTE NORMS
%==================================================
% Infinity norm
disp('=============== LINF NORM =================')
normRT_Linf = max(abs(inputRT(:)))
normKW_Linf = max(abs(sensor_data(:)))
quotient_Linf = normRT_Linf / normKW_Linf

% L1 norm
disp('=============== L1 NORM =================')
normRT_L1 = sum(abs(inputRT(:)))/length(inputRT(:))
normKW_L1 = sum(abs(sensor_data(:)))/length(sensor_data(:))
quotient_L1 = normRT_L1 / normKW_L1

% L2 norm
disp('=============== L2 NORM =================')
normRT_L2 = sqrt(sum((inputRT(:).*inputRT(:))/length(inputRT(:))))
normKW_L2 = sqrt(sum((sensor_data(:).*sensor_data(:))/length(sensor_data(:))))
quotient_L2 = normRT_L2 / normKW_L2




%============================================================================================================================================
% ADJOINT PROBLEM
%============================================================================================================================================

%==================================================
% kWave
%==================================================
p0_recon_PML = h5read('./output_data/Example69_adjoint_output_100sensors_homo.h5', '/p_final');
PML_size = 10;
pixelPressure_KW = max(0, p0_recon_PML(1+PML_size:end-PML_size, 1+PML_size:end-PML_size, 1+PML_size:end-PML_size));
plot_projection(pixelPressure_KW, 1);

%==================================================
% RT
%==================================================
% Import data
pixelPressureMatrix = importdata('./output_data/PixelPressure.dat', ' ', 0);
pixelPressure_RT = max(0, matrix2cube(pixelPressureMatrix, Nz));
plot_projection(pixelPressure_RT, 1);

%==================================================
% COMPUTE NORMS
%==================================================
% Infinity norm
disp('=============== LINF NORM =================')
normRT_Linf = max(abs(pixelPressure_RT(:)))
normKW_Linf = max(abs(pixelPressure_KW(:)))
quotient_Linf = normRT_Linf / normKW_Linf

% L1 norm
disp('=============== L1 NORM =================')
normRT_L1 = sum(abs(pixelPressure_RT(:)))/length(pixelPressure_RT(:))
normKW_L1 = sum(abs(pixelPressure_KW(:)))/length(pixelPressure_KW(:))
quotient_L1 = normRT_L1 / normKW_L1

% L2 norm
disp('=============== L2 NORM =================')
normRT_L2 = sqrt(sum((pixelPressure_RT(:).*pixelPressure_RT(:))/length(pixelPressure_RT(:))))
normKW_L2 = sqrt(sum((pixelPressure_KW(:).*pixelPressure_KW(:))/length(pixelPressure_KW(:))))
quotient_L2 = normRT_L2 / normKW_L2


