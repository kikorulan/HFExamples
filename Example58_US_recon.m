% Compute the rays and amplitude at the given domain
cd /cs/research/medim/projects2/projects/frullan/Documents/HighFreqCode/Examples/Ex58_US;

close all;
%clear all;
load gridLMatrix;
load prop_time_RT_126;


%====================================
% DOMAIN
%====================================
% Build domain
c0 = 1500;
factor = 0.1;
% Build Peaks
frame = 2.5;
x = -frame:(2*frame)/(grid.Ny-1):frame;
y = -frame:(2*frame)/(grid.Nx-1):frame;
[X, Y] = meshgrid(x, y);
p = peaks(X, Y);
cMatrix = c0*(ones(grid.Nx, grid.Ny) + factor*p/max(p(:)));
figure, surf(cMatrix, 'EdgeColor', 'none'), view(2), colorbar();

%==============================================================================================================
% SIMULATION
%==============================================================================================================

% Operators
K = grid.LMatrix;
Kadj = K';

% Initialization
x_0 = 1./reshape(1500*ones(grid.Nx, grid.Ny), [], 1);
x_k = x_0;
y_0 = prop_time_vec;
y_k = y_0;

% Parameters
tau = 1e3;
sigma = 1;
para.maxIter = 30;
para.constraint = 'positivity';
lambda = max(x_0)/1e4;
nIter = 200;

% Compute iterations
for i = 1:nIter
    % Prox g
    x_mat = reshape(x_k - tau*Kadj*y_k, grid.Nx, grid.Ny);
    if (mod(i, 10) == 0)
        figure, surf(1./x_mat, 'EdgeColor', 'none'), view(2), colorbar(), title(['Iter', num2str(i)]);
    end
    x_mat_1 = conTVdenoising(x_mat, lambda, para);
    x_k_1 = reshape(x_mat_1, [], 1);
    % Prox f
    y_k_1 = (y_k + sigma*K*(2*x_k_1 - x_k) - sigma*y_0)/(1+sigma);
    % Update
    x_k = x_k_1;
    y_k = y_k_1;
    pause(0.2);
    disp(i)
end

%==============================================================================================================
% VISUALIZATION
%==============================================================================================================

positionYBar = [700 700 500 400];
% Initial sound speed
figure;
surf(1e3*grid.xAxis, 1e3*grid.yAxis, cMatrix, 'EdgeColor', 'none');
axis([0 1e3*grid.xAxis(end) 0 1e3*grid.yAxis(end)]);
view(2);
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
box on;
set(gcf, 'pos', positionYBar);
%saveas(gcf, 'Example58_C', 'png');
%saveas(gcf, 'Example58_C.fig');

% Initial sound speed
figure;
surf(1e3*grid.xAxis, 1e3*grid.yAxis, 1./x_mat_1, 'EdgeColor', 'none');
axis([0 1e3*grid.xAxis(end) 0 1e3*grid.yAxis(end)]);
view(2);
xlabel('x [mm]');
ylabel('y [mm]');
colorbar();
caxis([1379 1650]);
box on;
set(gcf, 'pos', positionYBar);
%saveas(gcf, 'Example58_C_recon', 'png');
%saveas(gcf, 'Example58_C_recon.fig');
