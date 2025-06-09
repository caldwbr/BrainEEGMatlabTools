%Dyn PCA Flow



%% Load cleaned, upsampled voltage values
interpolatedVoltages = csvread('interpolatedVoltagesUnknownSongPlus.csv');
interpolatedVoltages = interpolatedVoltages(609900:3486300,:);
dataArray = interpolatedVoltages;

%% Create complex Morlet wavelet at 60 Hz
Fs = 1000; 
L = size(dataArray,1);
tEEG = (0:L-1) / Fs;
tCMW = -2:1/Fs:2;
frequencies = 60.0;%23.0:0.1:60.0;%[2.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0];
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


%% Convolve interpolatedVoltage w complex Morlet wavelet obtain analytic signal and from there the amplitude of 60 Hz
numChannels = size(dataArray,2);

% Initialize output: Amplitude per channel
ampEnvelopeAll = zeros(nData, numChannels);

% Do convolution for each channel
for ch = 1:numChannels
    eegData = dataArray(:,ch);
    dataX = fft(eegData, nConv);
    waveletX = fft(morletWavelets, nConv, 2);  % 2nd dim is time
    waveletX = waveletX ./ max(abs(waveletX)); % normalize kernel
    convResult = ifft(dataX .* waveletX.');
    convResult = convResult(halfK+1 : halfK+nData);
    ampEnvelopeAll(:,ch) = abs(convResult);  % Amplitude envelope for each channel
end

%% Find Dyn PCA weights and scores for each millisecond
%tf style sliding dynamic PCA
%PARAMETERS
Fs = 1000;
windowMs = 2000; % 2 sec window
halfWin = windowMs / 2;
sigmaMs = 500; % Gaussian sigma ~0.5 sec
timePoints = size(ampEnvelopeAll,1);
nChannels = size(ampEnvelopeAll,2);
nPCs = 3; % We want PC1-PC3

% Pad data (mirror padding)
ampPadded = padarray(ampEnvelopeAll, [halfWin 0], 'replicate', 'both');
timePointsPadded = size(ampPadded,1);

% Precompute Gaussian window
tVec = (-halfWin+1):(halfWin);
gaussWin = exp(-0.5 * (tVec/sigmaMs).^2);
gaussWin = gaussWin(:);
gaussWin = gaussWin / sum(gaussWin);

% Initialize output
PC_weights = zeros(timePoints, nChannels, nPCs); % time x channels x PCs

%% SLIDING WINDOW PCA
for t = 1:timePoints
    % Extract window
    dataWin = ampPadded(t:t+windowMs-1, :);
    
    % Apply Gaussian weighting to each channel
    weightedData = dataWin .* gaussWin;
    
    % Avoid log(0)
    weightedData(weightedData < 0.0001) = 0.0001;
    ampLog = log10(weightedData);
    
    % Z-score normalize
    ampZ = zscore(ampLog);
    
    % PCA
    [coeff, ~, ~] = pca(ampZ);
    
    % Store PC1-PC3 loadings, normalized
    for pcIdx = 1:nPCs
        pcVec = coeff(:, pcIdx);
        pcVec_norm = pcVec / norm(pcVec);
        PC_weights(t,:,pcIdx) = pcVec_norm';
    end
end

%% Convert to 3D trajectory like before
% Project original data onto time-resolved PCA weights
score_dynamic = zeros(timePoints, nPCs);

for t = 1:timePoints
    % Use ampZ from global prep (optional, or recalc local ampZ here)
    dataVec = log10(ampEnvelopeAll(t,:));
    dataVecZ = (dataVec - mean(dataVec)) ./ std(dataVec); 
    
    % Project onto time-local PCA basis
    for pcIdx = 1:nPCs
        pcVec = squeeze(PC_weights(t,:,pcIdx))';
        score_dynamic(t,pcIdx) = dot(dataVecZ, pcVec);
    end
end


%% STATIC 3D PLOT of dynamic PCA

figure('Position',[100 100 1000 800],'Color','w');
hold on;

plot3(score_dynamic(:,1), score_dynamic(:,2), score_dynamic(:,3), 'm-', 'LineWidth', 1.0);

xlabel('PC 1','FontSize',14);
ylabel('PC 2','FontSize',14);
zlabel('PC 3','FontSize',14);
title('Sliding PCA Trajectory','FontSize',16);
grid on;
axis equal;
box on;
view([30 20]);





%% Movie of using not using extra dimension for timeflow changes

trailDurationSec = 0.3; % how many seconds of trailing tail to show
trailLength = trailDurationSec * Fs;
videoFrameRate = 50; % frames per second
stepSize = round(Fs / videoFrameRate);
nFrames = floor(timePoints / stepSize);

% Get axis limits for consistent view
buffer = 0.2;
sampleWindow = score_dynamic(20*Fs:30*Fs,:); % same as before
pcMin = min(sampleWindow);
pcMax = max(sampleWindow);
pcRange = pcMax - pcMin;
pcMin = pcMin - buffer * pcRange;
pcMax = pcMax + buffer * pcRange;

%% CREATE MOVIE

figure('Position',[100 100 1000 800],'Color','w');
filename = 'entireGNBSsongsSession.mp4';
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = videoFrameRate;
open(v);

for iFrame = 1:nFrames
    clf;
    hold on;
    
    % Current time index
    tNow = (iFrame-1)*stepSize + 1;
    idxStart = max(1, tNow - trailLength);
    idxTrail = idxStart:tNow;
    
    % Plot full trajectory faint in background
    %plot3(score_dynamic(:,1), score_dynamic(:,2), score_dynamic(:,3), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 0.1);
    
    % Plot trail
    plot3(score_dynamic(idxTrail,1), score_dynamic(idxTrail,2), score_dynamic(idxTrail,3), 'b-', 'LineWidth', 2);
    
    % Current point
    plot3(score_dynamic(tNow,1), score_dynamic(tNow,2), score_dynamic(tNow,3), 'ro', 'MarkerFaceColor','r','MarkerSize',8);
    
    xlabel('PC 1','FontSize',14);
    ylabel('PC 2','FontSize',14);
    zlabel('PC 3','FontSize',14);
    title(sprintf('Sliding PCA Frame %d', iFrame), 'FontSize',16);
    
    grid on;
    axis equal;
    xlim([-2 2]);
    ylim([-5 5]);
    zlim([-4 4]);
    view([30 20]);
    box on;
    % (inside your for-loop)

    % Calculate rotating azimuth
    az0 = 30;
    el0 = 20;
    rotationPerFrame = 0.05;
    az = az0 + rotationPerFrame * iFrame;
    view([az el0]);

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);


%% MOTHERLODE: Dyn PCA but allowing present to update the recent past visually
%MOTHERLODE!!!!!!!
% --- PREP STAGE ---
startIdx = 609900;
endIdx = 3486300;
nHistory = 200;
stepSize = 20;
videoFrameRate = 50;

PC_weights_slice = PC_weights(startIdx:endIdx,:,:);
score_dynamic_slice = score_dynamic(startIdx:endIdx,:);

nTime = size(score_dynamic_slice, 1);
nFrames = floor(nTime / stepSize);

ampEnvelope_slice = ampEnvelopeAll(startIdx:endIdx,:);  % assuming ampEnvelope is [time x 16]

% Precompute adjusted scores
score_dynamic_adjusted = nan(nTime, nHistory, 3);

fprintf('Precomputing adjusted scores...\n');
for t = nHistory+1:nTime
    w = squeeze(PC_weights_slice(t,:,:));  % [16 x 3] weights at now
    for h = 1:nHistory
        sig = ampEnvelope_slice(t - h, :);  % [1 x 16]
        score_dynamic_adjusted(t, h, :) = sig * w;  % [1 x 3]
    end
end

% --- VIDEO STAGE ---
filename = 'BankWalk_Motherload_v1.mp4';
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = videoFrameRate;
open(v);

figure('Position',[100 100 1000 900],'Color','w');

% Axis limits
buffer = 0.1;
all_scores = reshape(score_dynamic_adjusted, [], 3);
xrange = range(all_scores(:,1), 'omitnan');
yrange = range(all_scores(:,2), 'omitnan');
zrange = range(all_scores(:,3), 'omitnan');

xlim_vals = [min(all_scores(:,1),[],'omitnan')-buffer*xrange, max(all_scores(:,1),[],'omitnan')+buffer*xrange];
ylim_vals = [min(all_scores(:,2),[],'omitnan')-buffer*yrange, max(all_scores(:,2),[],'omitnan')+buffer*yrange];
zlim_vals = [min(all_scores(:,3),[],'omitnan')-buffer*zrange, max(all_scores(:,3),[],'omitnan')+buffer*zrange];

fprintf('Generating video...\n');
for iFrame = 1:nFrames
    clf; hold on;

    tNow = iFrame * stepSize;

    % Plot history trail
    trail = squeeze(score_dynamic_adjusted(tNow, :, :));  % [200 x 3]
    plot3(trail(:,1), trail(:,2), trail(:,3), 'b-', 'LineWidth', 1.5);

    % Plot current tracer
    scatter3(score_dynamic_slice(tNow,1), score_dynamic_slice(tNow,2), score_dynamic_slice(tNow,3), ...
             50, 'r', 'filled');

    xlabel('PC1','FontSize',14);
    ylabel('PC2','FontSize',14);
    zlabel('PC3','FontSize',14);
    title(sprintf('t = %.3f sec', (tNow + startIdx)/1000), 'FontSize',16);
    grid on;
    xlim(xlim_vals); ylim(ylim_vals); zlim(zlim_vals);
    view([30 + 0.15*iFrame, 20]);
    box on;

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
