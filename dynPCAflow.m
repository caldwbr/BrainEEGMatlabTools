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

% SLIDING WINDOW PCA
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

% Convert to 3D trajectory like before
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

% CREATE MOVIE

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


%% Can save weights and scores
% --- Extract PC weights ---
PC1noGNBS = squeeze(PC_weights(:,:,1));  % [2876401 x 16]
PC2noGNBS = squeeze(PC_weights(:,:,2));  % [2876401 x 16]
PC3noGNBS = squeeze(PC_weights(:,:,3));  % [2876401 x 16]

% --- Save each PC weight matrix ---
writematrix(PC1noGNBS, 'PC1noGNBS.csv');
fprintf('Saved: PC1noGNBS.csv\n');

writematrix(PC2noGNBS, 'PC2noGNBS.csv');
fprintf('Saved: PC2noGNBS.csv\n');

writematrix(PC3noGNBS, 'PC3noGNBS.csv');
fprintf('Saved: PC3noGNBS.csv\n');

% --- Save score_dynamic ---
score_dynamicNoGNBS = squeeze(score_dynamic);  % [287401 x 3], just in case
writematrix(score_dynamicNoGNBS, 'score_dynamicNoGNBS.csv');
fprintf('Saved: score_dynamicNoGNBS.csv\n');



%% MOTHERLODE: Dyn PCA but allowing present to update the recent past visually
%MOTHERLODE!!!!!!!
% --- MOTHERLODE --- RETRO-STAMP MODE ---
startIdx = 2317600;
endIdx = 2657780;
nHistory = 400;  % 400 ms history
stepSize = 20;
videoFrameRate = 50;

PC_weights_slice = PC_weights(startIdx:endIdx,:,:);               % [time x 16 x 3]
score_dynamic_slice = score_dynamic(startIdx:endIdx,:);           % [time x 3]
ampEnvelope_slice = ampEnvelopeAll(startIdx:endIdx,:);            % [time x 16]

nTime = size(score_dynamic_slice, 1);
nFrames = floor(nTime / stepSize);

% Precompute adjusted scores using cast-added weight deltas
score_dynamic_adjusted = nan(nTime, nHistory, 3);  % [time x 400 ms x 3 PCs]

fprintf('Precomputing adjusted scores (cast-add deltas)...\n');
for t = 1 + nHistory : nTime
    for h = 1:nHistory
        t_h = t - h;
        % Start with the original weights at t-h
        w_h = squeeze(PC_weights_slice(t_h,:,:));  % [16 x 3]

        % Apply cast-add of all deltas from t_h+1 to t
        for shift = t_h+1 : t
            delta_w = squeeze(PC_weights_slice(shift,:,:) - PC_weights_slice(shift-1,:,:));
            w_h = w_h + delta_w;
        end

        amp_h = ampEnvelope_slice(t_h, :);  % [1 x 16]
        score_dynamic_adjusted(t, h, :) = amp_h * w_h;  % [1 x 3]
    end
end





%% PCA weights change visualization
filename = 'PC_weights_bargraph_labeled.mp4';
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = 50;
open(v);

figure('Position', [100 100 1200 600], 'Color', 'w');

nTime = size(PC_weights_slice, 1);
nElectrodes = 16;
nPCs = 3;

electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'Pz', 'Fz', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'Cz', 'P3', 'P4'};

% Build x-axis labels: 'Fp1-PC1', 'Fp1-PC2', ..., 'P4-PC3'
xLabels = cell(nElectrodes * nPCs, 1);
for i = 1:nElectrodes
    for j = 1:nPCs
        idx = (i-1)*nPCs + j;
        xLabels{idx} = sprintf('%s-PC%d', electrodeLabels{i}, j);
    end
end
x = 1:numel(xLabels);

% Optional: color-code PC1 = red, PC2 = green, PC3 = blue
barColors = repmat([1 0 0; 0 1 0; 0 0 1], nElectrodes, 1);  % [48 x 3]

for t = 1:20:nTime
    clf;

    weights = reshape(PC_weights_slice(t,:,:), [nElectrodes * nPCs, 1]);  % [48 x 1]

    b = bar(x, weights, 'FaceColor', 'flat');
    b.CData = barColors;

    ylim([-2 2]);  % Adjust based on actual range
    xlim([0 numel(xLabels)+1]);
    xticks(x);
    xticklabels(xLabels);
    xtickangle(90);
    title(sprintf('Timepoint %d (%.2f sec)', t, t / 50), 'FontSize', 16);
    ylabel('PCA Weight');

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);


% --- VIDEO GENERATION ---
filename = 'BankWalk_Motherload_v2_retrostamp.mp4';
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = videoFrameRate;
open(v);

figure('Position',[100 100 1000 900],'Color','w');

% Compute axis limits
buffer = 0.1;
flat_scores = reshape(score_dynamic_adjusted, [], 3);
xrange = range(flat_scores(:,1), 'omitnan');
yrange = range(flat_scores(:,2), 'omitnan');
zrange = range(flat_scores(:,3), 'omitnan');

xlim_vals = [min(flat_scores(:,1),[],'omitnan')-buffer*xrange, max(flat_scores(:,1),[],'omitnan')+buffer*xrange];
ylim_vals = [min(flat_scores(:,2),[],'omitnan')-buffer*yrange, max(flat_scores(:,2),[],'omitnan')+buffer*yrange];
zlim_vals = [min(flat_scores(:,3),[],'omitnan')-buffer*zrange, max(flat_scores(:,3),[],'omitnan')+buffer*zrange];

fprintf('Generating video...\n');
for iFrame = 1:nFrames
    clf; hold on;

    tNow = iFrame * stepSize;

    % Plot history trail using retrospective reprint with current weights
    trail = squeeze(score_dynamic_adjusted(tNow, :, :));  % [400 x 3]
    plot3(trail(:,1), trail(:,2), trail(:,3), 'b-', 'LineWidth', 1.5);

    % Plot current tracer dot (real PCA score at that moment)
    scatter3(score_dynamic_slice(tNow,1), score_dynamic_slice(tNow,2), score_dynamic_slice(tNow,3), ...
             50, 'r', 'filled');

    xlabel('PC1','FontSize',14);
    ylabel('PC2','FontSize',14);
    zlabel('PC3','FontSize',14);
    title(sprintf('t = %.3f sec', (tNow / 1000) + 609.9), 'FontSize',16);
    grid on;
    xlim(xlim_vals); ylim(ylim_vals); zlim(zlim_vals);
    view([30 + 0.15*iFrame, 20]);
    box on;

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Done!\n');




% --- SETTINGS ---
filename = 'PC_weightsDOODLE.mp4';
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = 50;
open(v);

figure('Position', [100 100 1200 600], 'Color', 'w');

% --- VARS ---
nElectrodes = 16;
nPCs = 3;
frameWindow = 1000;  % 20 seconds at 50 fps

% Electrode labels
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'Pz', 'Fz', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'Cz', 'P3', 'P4'};

% Build x-axis labels: 'Fp1-PC1', 'Fp1-PC2', ..., 'P4-PC3'
xLabels = cell(nElectrodes * nPCs, 1);
for i = 1:nElectrodes
    for j = 1:nPCs
        idx = (i-1)*nPCs + j;
        xLabels{idx} = sprintf('%s-PC%d', electrodeLabels{i}, j);
    end
end
x = 1:numel(xLabels);

% Color-code PC1 = red, PC2 = green, PC3 = blue
barColors = repmat([1 0 0; 0 1 0; 0 0 1], nElectrodes, 1);  % [48 x 3]

% --- FRAME LOOP ---
for t = 1:frameWindow
    clf;

    weights = reshape(PC_weights_slice(t,:,:), [nElectrodes * nPCs, 1]);  % [48 x 1]

    b = bar(x, weights, 'FaceColor', 'flat');
    b.CData = barColors;

    ylim([-0.8 0.8]);  % Adjust if needed
    xlim([0 numel(xLabels)+1]);
    xticks(x);
    xticklabels(xLabels);
    xtickangle(90);
    title(sprintf('Time %.2f sec (frame %d)', t/50, t), 'FontSize', 16);
    ylabel('PCA Weight');
    grid on;

    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);




%% Make plot of tracer speed
% === INPUTS ===
Fs = 1000;  % Sampling rate in Hz
dt = 1 / Fs;

% Assuming score_dynamic is [T x 3]
velocity_linear = nan(size(score_dynamic,1)-1,1);  % One less because diff

% Compute Euclidean distance between successive points
for t = 2:size(score_dynamic,1)
    delta = score_dynamic(t,:) - score_dynamic(t-1,:);
    velocity_linear(t-1) = norm(delta) / dt;  % mm/sec or whatever units
end

% === PLOT ===
tVec = (1:length(velocity_linear)) / Fs;

figure;
plot(tVec, velocity_linear, 'b');
xlabel('Time (s)');
ylabel('Linear Velocity in PCA Space');
title('Instantaneous Tracer Speed');
grid on;
