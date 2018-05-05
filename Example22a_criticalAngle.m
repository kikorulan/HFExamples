
% Measure computational time
tic;
start_time = clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INITIAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('workspace');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources and their respective Rays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
% Sources
for n = 1:nSources
    % Rays
    for j = 1:nRays
        plot(grid.source(n).x(1, :, j), grid.source(n).x(2, :, j), 'Color', colourMapV(n));
    end
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANALYSIS OF CRITICAL ANGLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial ray
nIni = 190;
pIni = grid.source(3).p(:, 1, nIni);
angleIni = atan(pIni(2)./pIni(1)) + pi*heaviside(-pIni(1));
% End ray
nEnd = 220;
pEnd = grid.source(3).p(:, 1, nEnd); 
angleEnd = atan(pEnd(2)./pEnd(1)) + pi*heaviside(-pEnd(1));

% Recompute rays for source 3
nRays = nEnd - nIni + 1;
x0 = grid.source(3).x0;

grid.newSource(x0, angleIni, angleEnd, nRays, step, tauMax);
grid.computeSourceQ(4);
grid.amplitudeODE(4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Verify superposition of solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Colour definition
nColours = 8;
colourList = cool(nColours);

figure;
hold on;
axis([0 Nx*dx 0 Ny*dy]);
% Rays
for j = nIni:nEnd
    plot(grid.source(3).x(1, :, j), grid.source(3).x(2, :, j), 'Color', 'g');
end
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(grid.source(4).x(1, :, j), grid.source(4).x(2, :, j), 'Color', colourList(colourNum, :));
end
legend('Sensor 1', 'Sensor 2', 'Sensor 3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Amplitude and Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amplitude
figure;
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    semilogy(grid.source(4).phi(:, :, j), grid.source(4).amplitude(:, :, j), 'Color', colourList(colourNum, :));
    hold on;
end

% Q
figure;
for j = 1:nRays
    colourNum = floor(nColours*(j-1)/nRays) + 1;
    plot(grid.source(4).phi(:, :, j), grid.source(4).q(:, :, j), 'Color', colourList(colourNum, :));
    hold on;
end


end_time = clock;
% Measure computational time
disp(['  total computation time ' num2str(etime(end_time, start_time))]);
