cd /cs/research/medim/projects2/projects/frullan/Documents/MATLAB/HighFreq/Examples/Ex26_timeStep_imp;
close all;

%load gridKWave.mat;
%load gridRT_time.mat;

inputRT = real(grid.source(1).pixelAPropagation);
inputKWave = gridKWave.signal;

%==============================
% Dimensions
%==============================
[xDimRT yDimRT zDimRT] = size(inputRT);
[xDimKWave yDimKWave zDimKWave] = size(inputKWave);

%==============================
% Norm calculation 
%==============================
vKWave = inputKWave(:);
maxValKWave = max(vKWave(vKWave < inf));
vRT = inputRT(:);
maxValRT = max(vRT(vRT < inf));

% L2 norm
L2RT = sqrt(sum(inputRT(:).^2));
L2KWave = sqrt(sum(inputKWave(:).^2));

% Maximum per frame
maxFrameRT = max(max(inputRT, [], 1), [], 2);
maxFrameKWave = max(max(inputKWave, [], 1), [], 2);
% L inf norm
LinfRT = repmat(maxFrameRT, [xDimRT yDimRT 1]);
LinfKWave = repmat(maxFrameKWave, [xDimKWave yDimKWave 1]);
% Plot inf norm
figure;
semilogy(permute(maxFrameRT, [3 2 1]), '-r');
hold on;
semilogy(permute(maxFrameKWave, [3 2 1]), '-b');
legend('RT', 'kWave');


% Assign chosen norm
normRT = LinfRT;
normKWave = LinfKWave;

%==============================
% Adjust Scaling
%==============================
inputRT_norm = inputRT./normRT;
inputKWave_norm = inputKWave./normKWave;

%==============================
% Synchronise signals
%==============================
frameKWave = 300;
inputKWave_norm_Rep = repmat(inputKWave_norm(:, :, 300), [1 1 zDimKWave]); 
% Compute difference
difference = (inputRT_norm - inputKWave_norm_Rep).*100;

% Compute synchronisation
dif2 = difference.^2;
sumDif = sum(sum(dif2));
syncError = frameKWave - find(sumDif == min(sumDif));
msg = strcat({'Synchronisation error KWave - RT: '}, int2str(syncError), {' frames'});
disp(msg{1});
if syncError > 0
    inputRT = inputRT(:, :, 1:end-syncError);
    inputKWave = inputKWave(:, :, syncError+1:end);
elseif syncError < 0
    inputRT = inputRT(:, :, -syncError+1:end);
    inputKWave = inputKWave(:, :, 1:end+syncError);
end

%%  %==============================
%%  %% Readjust scaling
%%  %==============================
%%  quotient = inputRT*normalise./inputKWave;
%%  
%%  scrollView(permute(fliplr(quotient), [2 1 3]), 3, [-5 5]);
%%  title('Quotient');
%% 
 
%==============================
% Error
%==============================
% Compute difference
difference = (inputRT_norm - inputKWave_norm).*100;

% Scroll View Error
scrollView(permute(fliplr(difference), [2 1 3]), 3, [-20 20]);
title('Error');

% Scroll View Log Error
logDifference = max(1e-5, abs(difference)/100);
scrollView(permute(fliplr(log10(logDifference)), [2 1 3]), 3, [-2 0]);
title('Log Error');

%==============================
% View Signals
%==============================
scrollView(permute(fliplr(inputRT_norm), [2 1 3]), 3, [-1 1]);
title('Ray Tracing');
scrollView(permute(fliplr(inputKWave_norm), [2 1 3]), 3, [-1 1]);
title('kWave');


% Difference with respect to 0
%scrollView(permute(fliplr(sensorReshape*100./factorKWave), [2 1 3]), 3, [-100 100]);

