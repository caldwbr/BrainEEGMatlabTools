% mindImage.m
% By Brad Caldwell, June 15, 2025.
% To take a mind image, 1. EEG, 2. CMW, 3. PCA.
% Need millisecond precision of each recording/processing step.
% It's a triple pass 'time/signal microscope' to convert brain activity to spatial schemas.

%% 1. Ms-resolution EEG. 
% 16-channel+ recommended. Load cleaned, upsampled voltage values. Pick a song to clip.
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'Pz', 'Fz', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'Cz', 'P3', 'P4'}; %Change if needed.
interpolatedVoltages = csvread('interpolatedVoltagesUnknownSongPlus.csv');
interpolatedVoltages = interpolatedVoltages(1484020:1672020,:);%Espresso
interpolatedVoltages = interpolatedVoltages(2927500:3267680,:);%Ordinary
interpolatedVoltages = interpolatedVoltages(893720:1212060,:);%TeenSpirit
interpolatedVoltages = interpolatedVoltages(3267680:end,:);%HeartAlone
dataArray = interpolatedVoltages;

%% 2. Ms-resolution CMW. 
% Create complex Morlet wavelet at frequenciers and then convolve interpolatedVoltage w complex Morlet wavelet obtain analytic signal and from there the amp or power of frequencies
Fs = 1000; 
L = size(dataArray,1);
tEEG = (0:L-1) / Fs;
tCMW = -5:1/Fs:5;
frequencies = 1.0:0.2:10.0;%[2.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0];
morletWavelets = zeros(length(frequencies), length(tCMW));
gaussianWindows = zeros(length(frequencies), length(tCMW));
sineWaves = zeros(length(frequencies), length(tCMW));
for i = 1:length(frequencies)
    f = frequencies(i);
    % Create complex sine wave
    sineWave = exp(1i * 2 * pi * f * tCMW);
    sineWaves(i, :) = sineWave; % Store the sine wave
    % Compute standard deviation for Gaussian based on 5 cycles
    s = 5 / (2 * pi * f);
    % Create Gaussian window
    gauss = exp(-tCMW.^2 / (2 * s^2));
    % Store Gaussian window
    gaussianWindows(i, :) = gauss;
    % Create complex Morlet wavelet
    morletWavelets(i, :) = sineWave .* gauss;
end

nData = length(tEEG);
nKern = length(morletWavelets(1,:));
nConv = nData + nKern - 1;
halfK = floor(nKern/2);

%Convolution bit
numChannels = size(dataArray,2);

ampEnvelopeAll = zeros(nData, numChannels, length(frequencies));  % [time x channels x freqs]

waveletX = fft(morletWavelets, nConv, 2);  % 2nd dim is time
waveletX = waveletX ./ max(abs(waveletX), [], 2);  % Normalize each wavelet

for fi = 1:length(frequencies)
    waveletFi = squeeze(waveletX(fi, :));
    
    for ch = 1:numChannels
        eegData = dataArray(:,ch);
        dataX = fft(eegData, nConv);
        convResult = ifft(dataX .* waveletFi.');
        convResult = convResult(halfK+1 : halfK+nData);
        ampEnvelopeAll(:,ch,fi) = abs(convResult);  % Store per freq
    end
end


%% 3. Ms-resolution PCA.
% Dynamic, millisecond-resolution PCA weights and scores
% No more interpolation!
% Dual guassians PCA good
timePoints = nData;
numFreqs = length(frequencies);
nPCs = 3;  % Only first 3 PCs
nCoarse = 224633;
coarseIdx = round(linspace(1, timePoints, nCoarse));
% Multi-Scale Gaussian Kernel (5 layers)
halfWin = round(3600 / 2);  % Centered window of ±900 ms
t = (-halfWin+1):(halfWin);

% Define sigmas
sigmas = [301, 31, 3];  % In ms
nGauss = numel(sigmas);
gaussStack = zeros(length(t), nGauss);

% Compute individual Gaussians
for i = 1:nGauss
    g = exp(-0.5 * (t / sigmas(i)).^2);
    g = g / sum(g);  % Normalize area
    gaussStack(:, i) = g;
end

% Average them together
gaussWin = mean(gaussStack, 2);
gaussWin = gaussWin / sum(gaussWin);  % Final normalization (just to be safe)

% halfWin = round(900 / 2);
% sigma1 = 220;
% sigma2 = 50;
% t = (-halfWin+1):(halfWin);
% 
% g1 = exp(-0.5 * (t / sigma1).^2);
% g2 = exp(-0.5 * (t / sigma2).^2);
% g1 = g1 / sum(g1);
% g2 = g2 / sum(g2);
% dualGauss = 0.5 * g1 + 0.5 * g2;
% %gaussWin = exp(-0.5 * ((-halfWin+1):(halfWin))'.^2 / sigmaMs^2);
% gaussWin = dualGauss;

% === INITIALIZE SCORES ONLY ===
PC_scores_1 = zeros(timePoints, numFreqs);
PC_scores_2 = zeros(timePoints, numFreqs);
PC_scores_3 = zeros(timePoints, numFreqs);

% === MAIN LOOP ===
for fi = 1:numFreqs
    disp(['Processing frequency: ' num2str(frequencies(fi)) ' Hz']);
    ampThisFreq = ampEnvelopeAll(:,:,fi);
    PC_weights_coarse = zeros(nCoarse, numChannels, nPCs);

    % COARSE PCA
    for i = 1:nCoarse
        if mod(i, round(nCoarse/100)) == 0
            fprintf('... %d%% done with coarse PCA\n', round(100 * i / nCoarse));
        end

        t = coarseIdx(i);
        tStart = max(1, t - halfWin + 1);
        tEnd = min(timePoints, t + halfWin);
        idxActual = (tStart:tEnd)';
        gwActual = gaussWin((tStart - t + halfWin):(tEnd - t + halfWin));

        dataWin = ampThisFreq(idxActual, :).^2;
        weightedData = dataWin .* gwActual;
        weightedData(weightedData < 0.0001) = 0.0001;
        dataLog = log10(weightedData);
        dataZ = zscore(dataLog);

        [coeff, ~, ~] = pca(dataZ);
        for pcIdx = 1:nPCs
            pcVec = coeff(:, pcIdx);
            pcVec_norm = pcVec / norm(pcVec);
            PC_weights_coarse(i,:,pcIdx) = pcVec_norm';
        end
    end


    % Precompute logZ data
    dataLogAll = log10(max(ampThisFreq, 0.0001));
    dataZAll = zscore(dataLogAll, 0, 2);

    % INLINE SCORE PROJECTION (no interpolation needed, we're doing PCA at every ms)
    for t = 1:timePoints
        dataVecZ = dataZAll(t, :);
        for pcIdx = 1:nPCs
            % === REPLACE INTERPOLATION WITH DIRECT ASSIGNMENT ===
            pcVec = squeeze(PC_weights_coarse(t,:,pcIdx));  % Already computed at t
            pcVec = pcVec / norm(pcVec);  % Just to be safe
            score = dot(dataVecZ, pcVec);
    
            if pcIdx == 1
                PC_scores_1(t,fi) = score;
            elseif pcIdx == 2
                PC_scores_2(t,fi) = score;
            elseif pcIdx == 3
                PC_scores_3(t,fi) = score;
            end
        end
    end

end

% === EXPORT TO BASE WORKSPACE ===
assignin('base', 'PC1', PC_scores_1);
assignin('base', 'PC2', PC_scores_2);
assignin('base', 'PC3', PC_scores_3);

% Save the PCA scores if desired
writematrix(PC_scores_1, 'PC1AlonePspectrumLowX.csv');
writematrix(PC_scores_2, 'PC2AlonePspectrumLowX.csv');
writematrix(PC_scores_3, 'PC3AlonePspectrumLowX.csv');

% Plot the movie, each band gets a PC1/2/3 tracer.
%% SETTINGS
Fs = 1000;
trailDurationSec = 0.22;
trailLength = round(trailDurationSec * Fs);
videoFrameRate = 50;
stepSize = round(Fs / videoFrameRate);
nFrames = floor(timePoints / stepSize);
%frequencies = freqs_01;
nFreqs = length(frequencies);


% % 🎨 Full 60-color "Duran Duran" colormap
% colors1 = [linspace(0, 0.2, 20)', linspace(0, 0.1, 20)', linspace(0.5, 1, 20)'];       % Dark blue to electric blue
% colors2 = [linspace(0.2, 0.9, 20)', linspace(0.1, 0.1, 20)', linspace(1, 0.6, 20)'];   % Blue to purplish pink
% colors3 = [linspace(0.9, 1, 20)', linspace(0.1, 0, 20)', linspace(0.6, 0.3, 20)'];     % Pink to hot magenta
% duranColormap = [colors1; colors2; colors3];  % 60x3 colormap

% % Teen Spirit 5-color grunge colormap (RGB 0–1)
% colorMap = [
%     0.102, 0.089, 0.071;  % shadow brown (lamp background)
%     0.600, 0.490, 0.298;  % golden tan (face highlight)
%     0.847, 0.712, 0.349;  % light yellow (sweater stripe)
%     0.259, 0.204, 0.133;  % earthy midtone (shirt & wall hue)
%     0.392, 0.157, 0.098   % tom drum red-brown (from this new image)
% ];
% colorIdx = mod(0:nFreqs-1, size(colorMap,1)) + 1;
% colorMap = colorMap(colorIdx, :);

% colorMap = [
%     0.113, 0.490, 0.741;  % vibrant pool blue
%     0.792, 0.792, 0.792;  % off-white newspaper
%     0.341, 0.717, 0.847;  % aqua cyan from baby highlight
%     0.019, 0.215, 0.396;  % deep pool shadow blue
%     0.717, 0.792, 0.678;  % greenish paper tint
%     0.866, 0.803, 0.360   % yellow from dollar bill
% ];
% colorMap = interp1(1:size(colorMap,1), colorMap, linspace(1, size(colorMap,1), 21), 'pchip');

% Make sure freqs_01 is defined and gives you 101 frequencies
%frequencies = freqs_01;
nFreqs = length(frequencies);  % Should be 101

% % Your base 5-color palette from image
% baseMap = [
%     0.004, 0.627, 0.792;  % vibrant pool blue
%     0.043, 0.345, 0.607;  % deep pool shadow blue
%     0.047, 0.424, 0.533;  % greenish paper tint
%     0.812, 0.745, 0.129;  % yellow from dollar bill
%     0.870, 0.850, 0.810   % offwhite newspaper
% ];
% 
% % Interpolate up to 101 colors
% colorMap = interp1(1:size(baseMap,1), baseMap, linspace(1, size(baseMap,1), nFreqs), 'pchip');
% 
% % Defensive check
% if size(colorMap,1) ~= nFreqs || size(colorMap,2) ~= 3
%     error('ColorMap interpolation failed: size is %dx%d, expected %dx3', ...
%           size(colorMap,1), size(colorMap,2), nFreqs);
% end

% % 🎯 Resample 60-color colormap to match 11 frequencies
% colorIdx = round(linspace(1, 60, nFreqs));
% colorMap = duranColormap(colorIdx, :);

% % === NEBULA HARDCODED ===
% nebula256 = [ % RGB rows, 0–1 fits carpenter fine
%     0.99, 0.01, 0.34
%     0.6, 0.16, 0.69
%     0.01, 0.53, 0.98
% ]; % You can extend this list to 256 rows (full nebula)

% nebula256 = [
%     0.22, 0.27, 0.42
%     0.85, 0.43, 0.33
%     0.95, 0.83, 0.71
% ];

nebula256 = [
    0,         0,    0.5000   % dark blue
    0,         0,    1.0000   % blue
    0,    1.0000,    1.0000   % cyan
    1.0000,    1.0000,    0   % yellow
    1.0000,         0,         0   % red

];

% Interpolate to match number of frequencies
colorMap = interp1(1:size(nebula256,1), nebula256, linspace(1, size(nebula256,1), nFreqs), 'pchip');

% Check
if size(colorMap,1) ~= nFreqs || size(colorMap,2) ~= 3
    error('ColorMap interpolation failed: got %dx%d, expected %dx3', ...
          size(colorMap,1), size(colorMap,2), nFreqs);
end


shiftAmount = 6;

% Load PC scores
PC1 = evalin('base', 'PC_scores_1');
PC2 = evalin('base', 'PC_scores_2');
PC3 = evalin('base', 'PC_scores_3');

%% Reorganize into scoresByFreq{f} = [PC1(t,f), PC2(t,f), PC3(t,f)]
scoresByFreq = cell(1, nFreqs);
for idx = 1:nFreqs
    f = frequencies(idx);
    scoresByFreq{idx} = [PC1(:,idx), PC2(:,idx), PC3(:,idx)];
end

% Offsets for spatial and time staggering
pc1Offsets = (10.0 - frequencies) * 0.15;
timeOffsets = (10.0 - frequencies) * 4;

%% Determine plot limits (log-transformed)
allTrail = [];
for idx = 1:nFreqs
    s = scoresByFreq{idx};
    s(:,1) = s(:,1) + pc1Offsets(idx);
    s = log10(s + shiftAmount);
    allTrail = [allTrail; s];
end
xMin = min(allTrail(:,1)); xMax = max(allTrail(:,1));
yMin = min(allTrail(:,2)); yMax = max(allTrail(:,2));
zMin = min(allTrail(:,3)); zMax = max(allTrail(:,3));

%% INIT FIGURE & VIDEO
figure('Color','k','Position',[100 100 1200 900]);
v = VideoWriter('heartAloneRingSpectralLowX.mp4', 'MPEG-4');
v.FrameRate = videoFrameRate;
open(v);



%% FRAME LOOP
for iFrame = 1:nFrames
    clf; hold on;
    tNow = (iFrame - 1) * stepSize + 1;

    for idx = 1:nFreqs
        scores = scoresByFreq{idx};
        color = colorMap(idx,:);
        tOffset = round(timeOffsets(idx));
        xOffset = pc1Offsets(idx);

        tUse = tNow - tOffset;
        if tUse < trailLength + 1 || tUse > size(scores,1)
            continue;
        end

        trailIdx = (tUse - trailLength + 1):tUse;
        trail = scores(trailIdx,:);
        trail(:,1) = trail(:,1) + xOffset;

        trailLog = log10(trail + shiftAmount);
        pt = scores(tUse,:); pt(1) = pt(1) + xOffset;
        ptLog = log10(pt + shiftAmount);

        % try
        %     plot3(trailLog(:,1), trailLog(:,2), trailLog(:,3), '-', ...
        %           'Color', color, 'LineWidth', 1.4);
        % catch ME
        %     disp('====== DEBUG INFO ======');
        %     disp(['size(trailLog) = ', mat2str(size(trailLog))]);
        %     disp(['size(scores) = ', mat2str(size(scores))]);
        %     disp(['idx = ', num2str(idx)]);
        %     disp(['color = ', mat2str(color)]);
        %     rethrow(ME)
        % end

        plot3(trailLog(:,1), trailLog(:,2), trailLog(:,3), '-', ...
              'Color', color, 'LineWidth', 1.4);
        plot3(ptLog(1), ptLog(2), ptLog(3), 'o', ...
              'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'MarkerSize', 6);
    end

    axis([xMin xMax yMin yMax zMin zMax]);
    set(gca, 'Color', 'k', 'Visible', 'off');
    view(30 + 0.25 * iFrame, 20);  % Slow z-axis spin

    writeVideo(v, getframe(gcf));
end

close(v);


