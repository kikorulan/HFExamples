%%% BUILD INITIAL PRESSURE

% Read data from files
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex86_3D_veins_het;
%cd /scratch0/NOT_BACKED_UP/frullan/Examples/Ex86_3D_veins_het;

%clear all;
close all;

% Read files
delimiterIn = ' ';
headerlinesIn = 0;

%========================================================================================================================
% DIMENSIONS
%========================================================================================================================
Nx = 240;
Ny = 240;
Nz = 240;

Nx_out = 240;
Ny_out = 240;
Nz_out = 240;
%========================================================================================================================
% BUILD PHANTOM
%========================================================================================================================
load veins/points_veins;
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
u0_p = zeros(Nx, Ny, Nz);
u0_p(pp_selec_lin) = 1;

disp('Veins computed');
% Convolve and threshold
uConv = u0_p;
para.maxIter = 4;
kernelSize = 2;
K = ones(kernelSize, kernelSize, kernelSize)/(kernelSize.^3);
uConv = convn(uConv, K, 'same');
thr1 = 0.3;
uConv(uConv <= thr1) = 0;
uConv(uConv >  thr1) = 1;
uConv = conTVdenoising(uConv, 0.1, para);
%thr2 = 0.5;
%uConv(uConv >  thr2) = 1;

% Save data
u0_resize = imresizen(uConv, Nx_out/Nx);
u0 = permute(u0_resize(:, :, floor(Nx_out/3) + 1:floor(2*Nx_out/3)), [3, 2, 1]);
plot_projection(u0, 1e-3);
%save input_data/initial_pressure.mat u0;
u0_matrix = cube2matrix(u0);
dlmwrite('./input_data/initial_pressure_veins_80x240x240.dat', u0_matrix, 'delimiter', ' ');

%========================================================================================================================
% BUILD SOUND SPEED
%========================================================================================================================
rng(1);
Nx_ss = floor(Nx/3)
gridR = gridRT(Nx_ss, 1, Ny, 1, Nz, 1);
% Set sound speed
c0 = 1580.0;
c = c0*ones(Nx_ss, Ny, Nz);
gridR.setCMatrix(c);
inc = 0.1;
gridR.randomiseC(.15, inc);

%==============================
% Plot sound speed
%==============================
plot_cube(gridR.c);
% Save data
c0_matrix = cube2matrix(gridR.c);
dlmwrite('./input_data/sound_speed.dat', c0_matrix, 'delimiter', ' ');

