close all;

load gridKWave_TS.mat;
%%  load gridKWave_IV.mat;
load grid_TS.mat;
%%  load grid_IV.mat;

inputRT = real(grid_TS.pixelAmplitudeSignal);
inputKWave = gridKWave_TS;

%==============================
%% Adjust scaling
%==============================
meanRT = mean(abs(inputRT(:)));
meanKWave = mean(abs(inputKWave(:)));
meanKWaveP = mean(abs(inputKWave(abs(inputKWave) > 0)));
normalise = meanKWave/meanRT;

%==============================
%% Synchronise signals
%==============================
frameKWave = 1200;
inputKWave_Rep = repmat(inputKWave(:, :, 1200), [1 1 size(inputKWave, 3)]); 
% Compute difference
difference = (inputRT*normalise - inputKWave_Rep)./...
                max(meanKWaveP, abs(inputKWave))*100;

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

%==============================
%% Readjust scaling
%==============================
quotient = inputRT*normalise./inputKWave;

scrollView(permute(fliplr(quotient), [2 1 3]), 3, [-5 5]);
title('Quotient');

%==============================
%% Error
%==============================
% Compute difference
difference = (real(inputRT)*normalise - inputKWave)./...
                max(meanKWaveP, abs(inputKWave))*100;

% Scroll View Error
scrollView(permute(fliplr(difference), [2 1 3]), 3, [-100 100]);
title('Error');

% Scroll View Log Error
logDifference = max(1e-5, abs(difference)/100);
scrollView(permute(fliplr(log10(logDifference)), [2 1 3]), 3, [-2 0]);
title('Log Error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% View Signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrollView(permute(fliplr(meanKWave*inputRT), [2 1 3]), 3, [-10*meanKWave*meanRT 10*meanKWave*meanRT]);
title('Ray Tracing');
scrollView(permute(fliplr(meanRT*inputKWave), [2 1 3]), 3, [-10*meanKWave*meanRT 10*meanKWave*meanRT]);
title('kWave');


% Difference with respect to 0
%scrollView(permute(fliplr(sensorReshape*100./meanKWave), [2 1 3]), 3, [-100 100]);

