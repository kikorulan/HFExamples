%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;

close all;
% Measure computational time
tic;
start_time = clock;

load recon_data_sensors.mat;
%==================================================================================
% SIMULATION
%==================================================================================
for m = 1:1
    msg = strcat({'Iteration '}, int2str(m), {'...'});
    disp(msg{1});
    clear pixelAReverse;
    parfor n = 1:256
        msg = strcat({'Source '}, int2str(n), {'...'});
        disp(msg{1});
        pixelAReverse{n} = grid.iterative_reconstruction(source(n));
        grid.pixelAReverse = grid.pixelAReverse + pixelAReverse{n};
    end

    for n = 1:256
        grid.pixelAReverse = grid.pixelAReverse + pixelAReverse{n};
    end
end
%==================================================================================
% VISUALISATION
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex31_reconstruction;
axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 400 600];
set(0,'DefaultFigurePaperPositionMode','auto');


%==============================
% Plot pressure
%==============================
pixelRT = max(0, grid.pixelAReverse)./max(grid.pixelAReverse(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelRT', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);

%========================================
% Reconstruction - kWave 
%========================================
%pixelKWave = p0_recon{257}/max(abs(p0_recon{257}(:)));
pixelKWave_TR = max(0, p0_recon{257});
pixelKWave_TR = pixelKWave_TR/max(pixelKWave_TR(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelKWave_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction - kWave (TR)');
saveas(gcf, 'Example31_kWave_recon_TR', 'png'); 
saveas(gcf, 'Example31_kWave_recon_TR.fig');

%========================================
% Error - RT
%========================================
errorRT_TR = pixelKWave_TR - pixelRT;
figure;
surf(X, Y, errorRT_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Error - Ray Tracing (TR)');
saveas(gcf, 'Example31_RT_error_TR', 'png'); 
saveas(gcf, 'Example31_RT_error_TR.fig');

%==============================
% Measure time
%==============================
end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);

% Save results
%save gridRT.mat grid nRays nSources x -v7.3;

cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples;

