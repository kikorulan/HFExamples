%================================================================================
% Example for gridRT class
%===============================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;

close all;
% Measure computational time
tic;
start_time = clock;

load recon_data.mat;
%==================================================================================
% SIMULATION
%==================================================================================
% Sequential execution
for m = 1:10
    msg = strcat({'Iteration '}, int2str(m), {'...'});
    disp(msg{1});
    for n = 764:-1:1
        msg = strcat({'Source '}, int2str(n), {'...'});
        disp(msg{1});
        grid.iterative_reconstruction(source(n));
    end
end


%%   Parallel execution
%%  or m = 1:3
%%     msg = strcat({'Iteration '}, int2str(m), {'...'});
%%     disp(msg{1});
%%     clear pixelAReverse;
%%     parfor n = 1:256
%%         msg = strcat({'Source '}, int2str(n), {'...'});
%%         disp(msg{1});
%%         pixelAReverse{n} = grid.iterative_reconstruction(source(n));
%%     end
%%  
%%     for n = 1:256
%%         grid.pixelAReverse = grid.pixelAReverse + pixelAReverse{n};
%%     end
%%  nd

%==================================================================================
% VISUALISATION
%==================================================================================
cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;
%cd /home/kiko/Documents/MATLAB/HighFreq/Examples/Ex36_reconstruction;
axisGrid = [0 (grid.Nx-1)*grid.dx 0 (grid.Ny-1)*grid.dy];

position = [700 700 400 600];
set(0,'DefaultFigurePaperPositionMode','auto');

%==============================
% Reconstruction - RT
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
title('Reconstruction RT - 10 iterations');
set(gcf, 'pos', position);
saveas(gcf, 'Example36_RT_recon_10iter', 'png');
saveas(gcf, 'Example36_RT_recon_10iter.fig');

%========================================
% Reconstruction - kWave
%========================================
pixelKWave = max(0, p0_recon)/max(p0_recon(:));
figure;
surf(grid.xAxis, grid.yAxis, pixelKWave', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Reconstruction k-Wave (TR)');
saveas(gcf, 'Example36_kWave_recon_TR', 'png');
saveas(gcf, 'Example36_kWave_recon_TR.fig');


%========================================
% Error - RT
%========================================
errorRT_TR = pixelKWave - pixelRT;
figure;
surf(grid.xAxis, grid.yAxis, errorRT_TR', 'EdgeColor', 'none');
view(2);
axis(axisGrid);
caxis([-.5, .5]);
colorbar();
box on;
xlabel('x (m)');
ylabel('y (m)');
set(gcf, 'pos', position);
title('Error RT - kWave (TR)');
saveas(gcf, 'Example36_RT_error_TR', 'png'); 
saveas(gcf, 'Example36_RT_error_TR.fig');

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

