%%% BUILD INITIAL PRESSURE

% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex63_3D_veins;
%clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%========================================================================================================================
% DIMENSIONS
%========================================================================================================================
Nx = 128;
Ny = 128;
Nz = 128;

%========================================================================================================================
% BUILD PHANTOM
%========================================================================================================================
load input_data/points_veins;
nPoints = size(PP, 1);
factor = 100;
inc = floor(nPoints/factor);

% Plot and select
pp_selection = [];
figure;
hold on;
for i = 1:factor
    index = i*inc + 1;
    index_1 = min((i+1)*inc, nPoints);
    if (i < 37  || (i > 43 && i < 47)||(i > 50 && i < 91))
        %scatter3(PP(index:index_1, 1), PP(index:index_1, 2), PP(index:index_1, 3));
        pp_selection = [pp_selection; PP(index:index_1, :)];
    end
end

%==============================
% Convert to Nx x Ny x Nz dimensions cube
%==============================
% Min to 0;
minX = min(pp_selection(:, 1));
pp_selection(:, 1) = pp_selection(:, 1) - minX;
minY = min(pp_selection(:, 2));
pp_selection(:, 2) = pp_selection(:, 2) - minY;
minZ = min(pp_selection(:, 3));
pp_selection(:, 3) = pp_selection(:, 3) - minZ;

% Center at cube of dim 1 x 1 x 1;
maxPP = max(pp_selection(:));
pp_selection = pp_selection/maxPP/1.5;
maxX = max(pp_selection(:, 1));
pp_selection(:, 1) = pp_selection(:, 1) + (1-maxX)/2;
maxY = max(pp_selection(:, 2));
pp_selection(:, 2) = pp_selection(:, 2) + (1-maxY)/2;
maxZ = max(pp_selection(:, 3));
pp_selection(:, 3) = pp_selection(:, 3) + (1-maxZ)/2;

% Plot
figure;
scatter3(pp_selection(:, 1), pp_selection(:, 2), pp_selection(:, 3));
view(2);

% Create grid
pp_selec_int(:, 1) = round(pp_selection(:, 1)*Nx);
pp_selec_int(:, 2) = round(pp_selection(:, 2)*Ny);
pp_selec_int(:, 3) = round(pp_selection(:, 3)*Nz);
pp_selec_lin = 1 + pp_selec_int(:, 1) + Nx*pp_selec_int(:, 2) + Nx*Ny*pp_selec_int(:, 3);
u0 = zeros(Nx, Ny, Nz);
u0(pp_selec_lin) = 1;


% Convolve and threshold
uConv = u0;
para.maxIter = 200;
kernelSize = 2;
K = ones(kernelSize, kernelSize, kernelSize)/(kernelSize.^3);
uConv = convn(uConv, K, 'same');
thr = 0.1;
uConv(uConv <= thr) = 0;
uConv(uConv >  thr) = 1;
uConv = conTVdenoising(uConv, 0.05, para);
uConv(uConv <= thr) = 0;
uConv(uConv >  thr) = 1;
plot_pixel(u0, 1, 1, false);
plot_pixel(uConv, 1, 1, false);

% Save data
u0 = uConv;
save input_data/initial_pressure.mat u0;
u0_matrix = cube2matrix(u0);
dlmwrite('input_data/initial_pressure_veins.dat', u0_matrix, 'delimiter', ' ');

%========================================================================================================================
% BUILD SOUND SPEED
%========================================================================================================================
%%  gridR = gridRT(Nx, 1, Ny, 1, Nz, 1);
%%  % Set sound speed
%%  c0 = 1500;
%%  c = c0*ones(Nx, Ny, Nz);
%%  gridR.setCMatrix(c);
%%  inc = 0.1;
%%  gridR.randomiseC(.15, inc);
%%  
%%  %==============================
%%  % Plot sound speed
%%  %==============================
%%  plot_cube(gridR.c);
%%  % Save data
%%  c0_matrix = cube2matrix(gridR.c);
%%  dlmwrite('input_data/sound_speed.dat', c0_matrix, 'delimiter', ' ');

