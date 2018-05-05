%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex56_RT_3D_4x4x4;

close all;
% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%==================================================
% Dimensions
%==================================================
% Import dimensions
dim = importdata('input_data/dimensions.dat', delimiterIn, headerlinesIn);
Nx = dim(1, 1);
dx = dim(2, 1);
Ny = dim(1, 2);
dy = dim(2, 2);
Nz = dim(1, 3);
dz = dim(2, 3);


%==================================================================================
% SIMULATION
%==================================================================================
% Import data
pixelPressureMatrix = importdata('output_data/PixelPressure.dat', delimiterIn, headerlinesIn);
pixelPressure = max(0, matrix2cube(pixelPressureMatrix, Nz));

% Define parameters
%%  lambda = 1e-1;
%%  lipschitz = 3;
%%  nIter = 5;
lambda = 1;
sigma = 1e4;
tau = 1e8;
theta = 1;
nEpochs = 1;
nSensors = 256;

% Parameters
para.maxIter = 20;
para.constraint = 'positivity';

% Proximal operators
proxG = @(x) conTVdenoising(x, 1e-2, para);
% Run iterations
pixelPressure_denoised = proxG(pixelPressure);


%==========================================================================================
% VISUALISATION
%==========================================================================================

%==================================================
% SCROLL RT
%==================================================
inputRT = pixelPressure;
inputRT(isinf(inputRT)) = 0;
sortedNormRT = sort(inputRT(:));
normRT = sortedNormRT(end-200);
inputRT_norm = inputRT/normRT;
scrollView(permute(fliplr(real(inputRT_norm)), [2 1 3]), 3, [0 1]);

%==================================================
% SCROLL RT - DENOISED
%==================================================
inputDE = pixelPressure_denoised;
inputDE(isinf(inputDE)) = 0;
sortedNormDE = sort(inputDE(:));
normDE = sortedNormDE(end-200);
inputDE_norm = inputDE/normDE;
scrollView(permute(fliplr(real(inputDE_norm)), [2 1 3]), 3, [0 1]);

%==================================================
% SCROLL RT - DENOISED - C++
%==================================================
% Import data
pixelPressureMatrix_cpp = importdata('output_data/image_denoised.dat', delimiterIn, 1);
pixelPressure_cpp = max(0, matrix2cube(pixelPressureMatrix_cpp.data, Nz));
inputRT_cpp = pixelPressure_cpp;
inputRT_cpp(isinf(inputRT_cpp)) = 0;
sortedNormRT_cpp = sort(inputRT_cpp(:));
normRT_cpp = sortedNormRT_cpp(end-200);
inputRT_norm_cpp = inputRT_cpp/normRT_cpp;
scrollView(permute(fliplr(real(inputRT_norm_cpp)), [2 1 3]), 3, [0 1]);

%==================================================
% SCROLL DIFFERENCE WITH C++
%==================================================
inputDif_norm = inputDE_norm-inputRT_norm_cpp;
scrollView(permute(fliplr(real(inputDif_norm)), [2 1 3]), 3, [-0.1 0.1]);

