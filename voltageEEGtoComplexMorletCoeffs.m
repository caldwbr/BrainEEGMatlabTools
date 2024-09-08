% Basic workflow to take voltage EEG and obtain an array
% of complex morlet wavelet coefficients for freqs 2/7/10/15/20/30/40hz
% With this array[freqs] of complex numbers, with just any
% timepoint, the complex number, you can obtain
% REAL real(eegConvress(freq,:,channel)
% IMAG imag(eegConvress(freq,:,channel)
% AMP abs(eegConvress(freq,:,channel)
% POWER abs(eegConvress(freq,:,channel).^2
% PHASE angle(eegConvress(freq,:,channel)
tableData = OpenBCIRAW20240509134543plusVapeNow;
tableData = tableData(6:end, :);
dataArray = table2array(tableData);
interpolatedVoltages = csvread('interpolatedVoltages.csv');
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4'};
dataArray = interpolatedVoltages;
Fs = 125; 
L = size(dataArray,1);
tEEG = (0:L-1) / Fs;
tCMW = -2:1/Fs:2;
frequencies = 5.0;%[2.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0];
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

eegDataX = fft(dataArray, nConv); 
eegKernXs = zeros(length(frequencies), nConv);
eegConvress = zeros(length(frequencies), nData, 16); 

for f = 1:length(frequencies)
    % FFT of the Morlet wavelet for the current frequency
    eegKernXs(f, :) = fft(morletWavelets(f, :), nConv);
    eegKernXs(f, :) = eegKernXs(f, :) ./ max(abs(eegKernXs(f, :))); % Normalize

    for ch = 1:16
        % Element-wise multiplication in frequency domain
        eegConvresX = eegDataX(:, ch) .* eegKernXs(f, :).'; % Ensure correct broadcasting
        % IFFT to get back to time domain
        eegConvresX = ifft(eegConvresX);
        % Truncate the result to remove edge artifacts
        eegConvress(f, :, ch) = eegConvresX(halfK + 1:end - halfK);
    end
end

%csvwrite('eegConvress.csv', eegConvress); %takes too long!

% Assuming eegConvress is a 7x449999x16 array
% 7 corresponds to the different frequency bands (2Hz, 7Hz, 10Hz, etc.)
% 449999 is the number of time points
% 16 corresponds to the channels

% Break eegConvress into individual frequency bands
complexNumbers2Hz = squeeze(eegConvress(1, :, :));  % Extract 2Hz data
complexNumbers5Hz = squeeze(eegConvress(1, :, :)); % Extract 5Hz data
complexNumbers7Hz = squeeze(eegConvress(2, :, :));  % Extract 7Hz data
complexNumbers10Hz = squeeze(eegConvress(3, :, :)); % Extract 10Hz data
complexNumbers15Hz = squeeze(eegConvress(4, :, :)); % Extract 15Hz data
complexNumbers20Hz = squeeze(eegConvress(5, :, :)); % Extract 20Hz data
complexNumbers30Hz = squeeze(eegConvress(6, :, :)); % Extract 30Hz data
complexNumbers40Hz = squeeze(eegConvress(7, :, :)); % Extract 40Hz data

csvwrite('complexNumbers2Hz.csv', complexNumbers2Hz);
csvwrite('complexNumbers5Hz.csv', complexNumbers5Hz);
csvwrite('complexNumbers7Hz.csv', complexNumbers7Hz);
csvwrite('complexNumbers10Hz.csv', complexNumbers10Hz);
csvwrite('complexNumbers15Hz.csv', complexNumbers15Hz);
csvwrite('complexNumbers20Hz.csv', complexNumbers20Hz);
csvwrite('complexNumbers30Hz.csv', complexNumbers30Hz);
csvwrite('complexNumbers40Hz.csv', complexNumbers40Hz);

% Split the complex data into real and imaginary parts
realParts = real(complexNumbers10Hz);    % Extract the real part
imaginaryParts = imag(complexNumbers10Hz); % Extract the imaginary part
% Save the real and imaginary parts as separate CSV files
csvwrite('complexNumbers10Hz_real.csv', realParts);
csvwrite('complexNumbers10Hz_imaginary.csv', imaginaryParts);


% Extract the first 5000 time points for channel 1 and channel 11
channel1Data = complexNumbers40Hz(1:5000, 1);
channel11Data = complexNumbers40Hz(1:5000, 11);

% Plot the real and imaginary parts of channel 1 vs channel 11
figure;

% Plot the real parts
subplot(2,1,1);
plot(real(channel1Data), 'b');
hold on;
plot(real(channel11Data), 'r');
title('Real Part of Channels 1 and 11 at 40Hz (First 5000 Time Points)');
xlabel('Time Points');
ylabel('Amplitude');
legend('Channel 1', 'Channel 11');
hold off;

% Plot the imaginary parts
subplot(2,1,2);
plot(imag(channel1Data), 'b');
hold on;
plot(imag(channel11Data), 'r');
title('Imaginary Part of Channels 1 and 11 at 40Hz (First 5000 Time Points)');
xlabel('Time Points');
ylabel('Amplitude');
legend('Channel 1', 'Channel 11');
hold off;




%Interactive application to show helix at 16 ch.
% Data details
samplingRate = 125; % 125 Hz sampling rate
windowSize = 5 * samplingRate; % 5 seconds window (625 time points)
totalTime = 3600; % 1 hour in seconds

% Create figure and set up key press detection
figh = figure('Position', [0, 0, 1400, 800], 'Color', 'w');
isPlaying = false; % Flag to indicate if playing or paused

% Set key press function
set(figh, 'KeyPressFcn', @(src, event) keyPress(event));

% Function to handle key press
function keyPress(event)
    if strcmp(event.Key, 'space')
        isPlaying = ~isPlaying; % Toggle play/pause
    end
end

% Loop through the data
t = 1; % Start time point
while t < totalTime * samplingRate - windowSize
    
    % If the user has pressed space to play
    if isPlaying
        clf; % Clear figure
        
        for ch = 1:16
            % Extract the data for the current channel over the window
            currentData = complexNumbers10Hz(t:(t + windowSize - 1), ch);

            % Plot helix in complex plane or sine wave (choose one)
            subplot(4, 4, ch); % Place channels in a 4x4 grid

            % Option 1: Helix in complex plane
            plot3(real(currentData), imag(currentData), (1:windowSize) / samplingRate, 'b');
            grid on;
            title(sprintf('Channel %d', ch));
            xlabel('Real');
            ylabel('Imaginary');
            zlabel('Time (s)');
            
            % Option 2: Real axis sine wave (comment out helix above if using sine wave)
            % plot((1:windowSize) / samplingRate, real(currentData), 'b');
            % title(sprintf('Channel %d', ch));
            % xlabel('Time (s)');
            % ylabel('Amplitude');
        end

        % Update the frame and wait to simulate real time
        drawnow;
        t = t + windowSize; % Move to the next time window
        
        % Simulate real-time playback (adjust sleep time as needed)
        pause(5 / samplingRate); % Adjust playback rate
    end
end




% Preallocate arrays for real, imaginary, and phase components
realParts = real(eegConvress);
imagParts = imag(eegConvress);
phaseParts = angle(eegConvress); % Precompute phase values
phaseDegParts = rad2deg(phaseParts);


% Save these arrays to avoid repeated calculations
save('realParts.mat', 'realParts', '-v7.3');
save('imagParts.mat', 'imagParts', '-v7.3');
save('phaseParts.mat', 'phaseParts', '-v7.3');

% Save as CSV files
csvwrite('realParts.csv', realParts);
csvwrite('imagParts.csv', imagParts);
csvwrite('phaseParts.csv', phaseParts);
csvwrite('phaseDegParts.csv', phaseDegParts);
csvwrite('phaseDiffDegParts.csv', phaseDiffDegParts);
csvwrite('interpolatedVoltages.csv', interpolatedVoltages);

% 2 Hz
phaseDiffDeg2Hz = squeeze(phaseDiffDegParts(1, :, :));
csvwrite('phaseDiffDeg2Hz.csv', phaseDiffDeg2Hz);

% 7 Hz
phaseDiffDeg7Hz = squeeze(phaseDiffDegParts(2, :, :));
csvwrite('phaseDiffDeg7Hz.csv', phaseDiffDeg7Hz);

% 10 Hz
phaseDiffDeg10Hz = squeeze(phaseDiffDegParts(3, :, :));
csvwrite('phaseDiffDeg10Hz.csv', phaseDiffDeg10Hz);

% 15 Hz
phaseDiffDeg15Hz = squeeze(phaseDiffDegParts(4, :, :));
csvwrite('phaseDiffDeg15Hz.csv', phaseDiffDeg15Hz);

% 20 Hz
phaseDiffDeg20Hz = squeeze(phaseDiffDegParts(5, :, :));
csvwrite('phaseDiffDeg20Hz.csv', phaseDiffDeg20Hz);

% 30 Hz
phaseDiffDeg30Hz = squeeze(phaseDiffDegParts(6, :, :));
csvwrite('phaseDiffDeg30Hz.csv', phaseDiffDeg30Hz);

% 40 Hz
phaseDiffDeg40Hz = squeeze(phaseDiffDegParts(7, :, :));
csvwrite('phaseDiffDeg40Hz.csv', phaseDiffDeg40Hz);

csvwrite('phaseDiffDeg2Hz.csv', phaseDiffDeg2Hz);
csvwrite('phaseDiffDeg7Hz.csv', phaseDiffDeg7Hz);
csvwrite('phaseDiffDeg10Hz.csv', phaseDiffDeg10Hz);
csvwrite('phaseDiffDeg15Hz.csv', phaseDiffDeg15Hz);
csvwrite('phaseDiffDeg20Hz.csv', phaseDiffDeg20Hz);
csvwrite('phaseDiffDeg30Hz.csv', phaseDiffDeg30Hz);
csvwrite('phaseDiffDeg40Hz.csv', phaseDiffDeg40Hz);



phaseDiffDeg2Hz = phaseDiffDeg2Hz(1:5:449999,:);
phaseDiffDeg7Hz = phaseDiffDeg7Hz(1:5:449999,:);
phaseDiffDeg10Hz = phaseDiffDeg10Hz(1:5:449999,:);
phaseDiffDeg15Hz = phaseDiffDeg15Hz(1:5:449999,:);
phaseDiffDeg20Hz = phaseDiffDeg20Hz(1:5:449999,:);
phaseDiffDeg30Hz = phaseDiffDeg30Hz(1:5:449999,:);
phaseDiffDeg40Hz = phaseDiffDeg40Hz(1:5:449999,:);

% Save comparison labels to a CSV file
writecell(comparisonLabels, 'comparisonLabels.csv');

% Define the filename for the .mat file
vapeComplexNumbers2to40Hz = 'eegConvress.mat';

% Save the variable eegConvress to the .mat file using version 7.3
save(vapeComplexNumbers2to40Hz, 'eegConvress', '-v7.3');

hz = linspace(0,srate/2,floor(nConv/2)+1);
figure(9), clf
subplot(211), hold on
plot(hz,abs(dataX(1:length(hz))),'b')
plot(hz,abs(kernX(1:length(hz))).*max(abs(dataX(:,1).'))/2)
plot(hz,abs(convresX(1:length(hz))),'k','linew',2)
set(gca,'xlim',[0 frex*2])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
legend({'Data spectrum';'Wavelet spectrum';'Convolution result'})

subplot(212), hold on
plot(t,centeredSamples(:,1),'b')
plot(t,real(convres),'k','linew',2)
%set(gca,'xlim',[-.1 1.3])
legend({'LFP data';'Convolution result'})
xlabel('Time (s)'), ylabel('Activity (\muV')





% Length of data
L = size(centeredSamples, 1);

% Time vector
t = (0:L-1) / Fs;

% Compute FFT of each channel
Y = fft(centeredSamples);

% Frequency axis
f = Fs*(0:(L/2))/L;

% Find index closest to 10 Hz
[~, idx_10Hz] = min(abs(f - 10));

% Create sine waves from FFT components
figure;
for channel = 1:size(centeredSamples, 2)
    % Get amplitude and phase at 10 Hz
    amplitude_10Hz = 2 * abs(Y(idx_10Hz+1, channel)) / L;  % Single-sided amplitude
    phase_10Hz = angle(Y(idx_10Hz+1, channel));  % Phase
    
    % Generate sine wave at 10 Hz with correct amplitude and phase
    sineWave = amplitude_10Hz * sin(2*pi*10*t + phase_10Hz);
    
    % Plot the sine wave for the current channel
    subplot(size(centeredSamples, 2), 1, channel);
    plot(t, sineWave, 'b');  % Plotting in blue color
    title(['Channel ' num2str(channel) ' - ' electrodeLabels{channel} ' - 10 Hz Sine Wave']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([t(1) t(end)]);  % Adjust xlim based on your data length
end

% Display the first 20 values for each of the sixteen channels
disp('First 20 values of centeredData for each channel:');
disp(centeredSamples(1:20, :));

% Extract the first 250 time points
sampleSection = centeredData(1000:5000, :);

meanSamples = mean(sampleSection, 1);

centeredSamples = sampleSection - meanSamples;

% Assuming centeredSamples contains your centered EEG data for the sample section

% Electrode labels corresponding to each channel
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4'};

% Sampling rate
Fs = 125;  % Adjust as per your actual sampling rate

% Length of data
L = size(centeredSamples, 1);

% Time vector
t = (0:L-1) / Fs;

% Compute FFT of each channel
Y = fft(centeredSamples);

% Frequency axis
f = Fs*(0:(L/2))/L;

% Find index closest to 10 Hz
[~, idx_10Hz] = min(abs(f - 10));

% Create sine waves from FFT components
figure;
for channel = 1:size(centeredSamples, 2)
    % Get amplitude and phase at 10 Hz
    amplitude_10Hz = 2 * abs(Y(idx_10Hz+1, channel)) / L;  % Single-sided amplitude
    phase_10Hz = angle(Y(idx_10Hz+1, channel));  % Phase
    
    % Generate sine wave at 10 Hz with correct amplitude and phase
    sineWave = amplitude_10Hz * sin(2*pi*10*t + phase_10Hz);
    
    % Plot the sine wave for the current channel
    subplot(size(centeredSamples, 2), 1, channel);
    plot(t, sineWave, 'b');  % Plotting in blue color
    title(['Channel ' num2str(channel) ' - ' electrodeLabels{channel} ' - 10 Hz Sine Wave']);
    xlabel('Time (s)');
    ylabel(electrodeLabels{channel});
    xlim([t(1) t(end)]);  % Adjust xlim based on your data length
end








%Complex morlet wavelet
srate = 125;
time = -1:1/srate:1;
frex = 10.0;
frex7 = 7.0;
frex2 = 2.0;
frex4 = 4.0;
frex40 = 40.0;
csw = exp( 1i*2*pi*frex*time );
csw2 = exp( 1i*2*pi*frex2*time );
csw4 = exp( 1i*2*pi*frex4*time );
csw7 = exp( 1i*2*pi*frex7*time );
csw40 = exp( 1i*2*pi*frex40*time );
fwhm = .5; %width of gaussian in seconds
gaus_win = exp( (-4*log(2)*time.^2) / fwhm^2 );
cmw = csw .* gaus_win;
cmw7 = csw7 .* gaus_win;
cmw2 = csw2 .* gaus_win;
cmw4 = csw4 .* gaus_win;
cmw40 = csw40 .* gaus_win;

figure(8), clf
subplot(211), hold on
plot(time,real(cmw),'b')
plot(time,imag(cmw),'r--')
xlabel('Time (s)'), ylabel('Amplitude')
legend({'real part';'imag part'})
title('Complex Morlet wavelet in the time domain')

subplot(212)
plot3(time, real(cmw), imag(cmw),'k','linew',3)
xlabel('Time (s)'), ylabel('Real'), zlabel('Imag')

nData = length(t);
nKern = length(cmw);
nConv = nData + nKern - 1;
halfK = floor(nKern/2);

dataX = fft(centeredSamples,nConv);
kernX = fft(cmw,nConv);
kernX7 = fft(cmw7,nConv);
kernX7 = kernX7 ./ max(kernX7);
kernX = kernX ./ max(kernX);
kernX2 = fft(cmw2,nConv);
kernX2 = kernX2 ./ max(kernX2);
kernX4 = fft(cmw4,nConv);
kernX4 = kernX4 ./ max(kernX4);
kernX40 = fft(cmw40,nConv);
kernX40 = kernX40 ./ max(kernX40);

convresX = dataX .* kernX.';
convresX7 = dataX .* kernX7.';
convres7 = ifft(convresX7);
convres7 = convres7(halfK+1:end-halfK,:);
convres = ifft(convresX);
convres = convres(halfK+1:end-halfK,:);

convresX2 = dataX .* kernX2.';
convres2 = ifft(convresX2);
convres2 = convres2(halfK+1:end-halfK,:);
convresX4 = dataX .* kernX4.';
convres4 = ifft(convresX4);
convres4 = convres4(halfK+1:end-halfK,:);
convresX40 = dataX .* kernX40.';
convres40 = ifft(convresX40);
convres40 = convres40(halfK+1:end-halfK,:);

hz = linspace(0,srate/2,floor(nConv/2)+1);
figure(9), clf
subplot(211), hold on
plot(hz,abs(dataX(1:length(hz))),'b')
plot(hz,abs(kernX(1:length(hz))).*max(abs(dataX(:,1).'))/2)
plot(hz,abs(convresX(1:length(hz))),'k','linew',2)
set(gca,'xlim',[0 frex*2])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
legend({'Data spectrum';'Wavelet spectrum';'Convolution result'})

subplot(212), hold on
plot(t,centeredSamples(:,1),'b')
plot(t,real(convres),'k','linew',2)
%set(gca,'xlim',[-.1 1.3])
legend({'LFP data';'Convolution result'})
xlabel('Time (s)'), ylabel('Activity (\muV')


% Define colors for each channel
colors = lines(16); % You can use any colormap you prefer

figure;
hold on;
for i = 1:16
    plot3(t, real(convres(:,i)), imag(convres(:,i)), 'Color', colors(i,:), 'LineWidth', 2);
end
hold off;

% Add labels for each channel
legend('Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', 'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4');
xlabel('Time');
ylabel('Real part');
zlabel('Imaginary part');
title('Complex Wavelet Coefficients for EEG Channels');

figure
plot3(t, real(convres2(:,1)), imag(convres2(:,1)),'r','linew',2), hold on
plot3(t, real(convres4(:,1)), imag(convres4(:,1)),'y','linew',2)
plot3(t, real(convres7(:,1)), imag(convres7(:,1)),'g','linew',2)
plot3(t, real(convres(:,1)), imag(convres(:,1)),'b','linew',2)
plot3(t, real(convres40(:,1)), imag(convres40(:,1)),'m','linew',2)
legend({'2 Hz at Fp1';'4 Hz at Fp1';'7 Hz at Fp1';'10 Hz at Fp1';'40 Hz at Fp1'})
xlabel('Time (s)'), ylabel('Real'), zlabel('Imag')
title('Complex Signal from Sedated State EEG')



% Parameters
eeg_signal = centeredSamples(:, 1);  % Use only the first channel
timeGPT = -2:1/srate:2;  % Time vector for wavelet
freqs = 2:0.5:62.5;  % Frequency range from 2 Hz to 62.5 Hz
num_cycles = 6;  % Number of cycles for the Morlet wavelet


% Initialize the combined signal as a complex array
combined_signal = complex(zeros(size(eeg_signal)));

% Compute the wavelet coefficients and sum them
for fi = 1:length(freqs)
    wavelet = morlet_wavelet(freqs(fi), srate, num_cycles, timeGPT);
    conv_result = conv(eeg_signal, wavelet, 'same');
    combined_signal = combined_signal + conv_result(1:length(eeg_signal));
end

% Visualize the combined complex signal as a helix
figure;
hold on;

real_part = real(combined_signal);
imag_part = imag(combined_signal);
time_vector = (0:length(eeg_signal)-1) / srate;

plot3(time_vector, real_part, imag_part, 'LineWidth', 2);

xlabel('Time (s)');
ylabel('Real Part');
zlabel('Imaginary Part');
title('Helix Visualization of Combined EEG Signal in Complex Plane');
grid on;
view(3);
hold off;


slotA = 2 +5i;
slotB = 1 - 2i;
slotC = -8 +2.3i;
slotD = -1 - 1.5i;
slotE = 3+2i;
slotF = 4+3i;
slotG = 0.2+3i;
slotH = -3 +1.5i;
slotI = -1.8 -2i;

% Create the 3x3 array
complexArray = [slotA, slotB, slotC;
                slotD, slotE, slotF;
                slotG, slotH, slotI];

% Save the complex array to a .mat file
save('complexArray.mat', 'complexArray');






% Number of channels
nChannels = length(electrodeLabels);

% Preallocate the phase difference matrix
[nFrequencies, nTimepoints, ~] = size(phaseDegParts);
nPairs = nchoosek(nChannels, 2); % Calculate number of unique pairs
phaseDiffDegParts = zeros(nFrequencies, nTimepoints, nPairs);
comparisonLabels = cell(nPairs, 1);

% Generate all unique pairs
pairs = nchoosek(1:nChannels, 2);

% Calculate phase differences for each pair and create labels
for i = 1:nPairs
    ch1 = pairs(i, 1);
    ch2 = pairs(i, 2);
    phaseDiffDegParts(:, :, i) = phaseDegParts(:, :, ch1) - phaseDegParts(:, :, ch2);
    comparisonLabels{i} = sprintf('%s:%s', electrodeLabels{ch1}, electrodeLabels{ch2});
end

%NECESSARY TO KEEP PHASE BETWEEN -180 and +180!! OTHERWISE IT GOES TO -360
%to +360!!!
phaseDiffDegParts(phaseDiffDegParts > 180) = phaseDiffDegParts(phaseDiffDegParts > 180) - 360;
phaseDiffDegParts(phaseDiffDegParts < -180) = phaseDiffDegParts(phaseDiffDegParts < -180) + 360;

% Display a few examples of the labels to verify
disp('Example comparison labels:');
disp(comparisonLabels(1:5));



% Given dataArray (448209 x 17 double array for large EEG data)
% Columns 1-16: Voltages, Column 17: Timepoints (unix)

% Extract voltages (448209 x 16) and timepoints (448209 x 1)
voltages = dataArray(:, 1:16);
timepoints = dataArray(:, 30);

% Convert timepoints to milliseconds
timepoints_ms = (timepoints - min(timepoints)) * 1000;

% Define the desired sampling rate for extrapolated values
desired_interval_ms = 8; % 8 ms interval

% Calculate the number of points for the interpolatedVoltages array
num_points = ceil((max(timepoints_ms) - min(timepoints_ms)) / desired_interval_ms) + 1;
interpolatedVoltages = zeros(num_points, 16);

% Define chunk size (number of points per chunk)
chunk_size = 100000; % Adjust based on your system's memory capacity

% Process the data in chunks for each channel
for channel = 1:16
    interpolatedChannel = [];
    for start_idx = 1:chunk_size:length(timepoints_ms)
        end_idx = min(start_idx + chunk_size - 1, length(timepoints_ms));
        
        % Get the current chunk
        chunk_timepoints = timepoints_ms(start_idx:end_idx);
        chunk_voltages = voltages(start_idx:end_idx, channel);
        
        % Ensure unique timepoints for the chunk
        [unique_timepoints, unique_indices] = unique(chunk_timepoints);
        unique_voltages = chunk_voltages(unique_indices);
        
        % Create the time vector for the smooth curve
        t_smooth = linspace(min(unique_timepoints), max(unique_timepoints), ...
                            ceil((max(unique_timepoints) - min(unique_timepoints)) / desired_interval_ms) + 1);
        
        % Fit a smooth curve using spline interpolation
        voltage_smooth = interp1(unique_timepoints, unique_voltages, t_smooth, 'spline');
        
        % Create the time vector for the extrapolated values within this chunk
        t_extrapolated = unique_timepoints(1):desired_interval_ms:unique_timepoints(end);
        
        % Query the smooth curve for the extrapolated values
        voltage_extrapolated = interp1(t_smooth, voltage_smooth, t_extrapolated, 'spline');
        
        % Append the extrapolated values to the result for this channel
        interpolatedChannel = [interpolatedChannel; voltage_extrapolated(:)];
    end
    
    % Store the result for this channel
    interpolatedVoltages(1:length(interpolatedChannel), channel) = interpolatedChannel;
end

% Ensure that interpolatedVoltages is trimmed to the correct number of points
interpolatedVoltages = interpolatedVoltages(1:num_points, :);

% Now interpolatedVoltages contains the extrapolated values for all channels



figure; plot(timepoints(1000:1100), voltages(1000:1100, 1), 'ro', timepoints(1000:1100), interpolatedVoltages(1000:1100, 1), 'b-');
xlabel('Time'); ylabel('Voltage'); legend('Original Data', 'Interpolated Data');


% Define the range of interest
range_start = 100;
range_end = 200100;

% Create the figure
figure;
hold on;

% Plot the original voltages
plot(timepoints_ms(range_start:range_end), voltages(range_start:range_end, 1), 'ro', 'MarkerSize', 8, 'DisplayName', 'Original Data');

% Plot the interpolated voltages as a blue line
% Match the indices of full_t_extrapolated to timepoints_ms for alignment
[~, idx_start] = min(abs(full_t_extrapolated - timepoints_ms(range_start)));
[~, idx_end] = min(abs(full_t_extrapolated - timepoints_ms(range_end)));

plot(full_t_extrapolated(idx_start:idx_end), interpolatedVoltages(idx_start:idx_end, 1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Interpolated Data');

% Plot the interpolated voltages as large black triangles
plot(full_t_extrapolated(idx_start:idx_end), interpolatedVoltages(idx_start:idx_end, 1), 'k^', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Interpolated Points');

% Add labels and legend
xlabel('Time (ms)');
ylabel('Voltage');
legend show;
title('Original and Interpolated Voltages');
grid on;
hold off;











%Phase diffs array July 12th attempt

%load
phaseDiffDeg2Hz = csvread('phaseDiffDeg2Hz25fps.csv');
phaseDiffDeg7Hz = csvread('phaseDiffDeg7Hz25fps.csv');
phaseDiffDeg10Hz = csvread('phaseDiffDeg10Hz25fps.csv');
phaseDiffDeg15Hz = csvread('phaseDiffDeg15Hz25fps.csv');
phaseDiffDeg20Hz = csvread('phaseDiffDeg20Hz25fps.csv');
phaseDiffDeg30Hz = csvread('phaseDiffDeg30Hz25fps.csv');
phaseDiffDeg40Hz = csvread('phaseDiffDeg40Hz25fps.csv');

% Define specific left/right and front/back comparisons
leftRightComparisons = {
    'Fp1', 'Fp2';
    'F7', 'F8';
    'F3', 'F4';
    'T3', 'T4';
    'C3', 'C4';
    'T5', 'T6';
    'P3', 'P4';
    'O1', 'O2'
};

frontBackComparisons = {
    'Fp1', 'O1';
    'Fp2', 'O2';
    'F7', 'T5';
    'F8', 'T6';
    'F3', 'P3';
    'F4', 'P4'
};

% Initialize data arrays
data = {phaseDiffDeg2Hz, phaseDiffDeg7Hz, phaseDiffDeg10Hz, ...
        phaseDiffDeg15Hz, phaseDiffDeg20Hz, phaseDiffDeg30Hz, ...
        phaseDiffDeg40Hz};

% Define channels and colors
channels = {'Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', ...
            'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4'};
channelPairs = {'Fp1-Fp2', 'Fp1-C3', 'Fp1-C4', 'Fp1-T5', 'Fp1-T6', ...
                'Fp1-O1', 'Fp1-O2', 'Fp1-F7', 'Fp1-F8', 'Fp1-F3', ...
                'Fp1-F4', 'Fp1-T3', 'Fp1-T4', 'Fp1-P3', 'Fp1-P4', ...
                'Fp2-C3', 'Fp2-C4', 'Fp2-T5', 'Fp2-T6', 'Fp2-O1', ...
                'Fp2-O2', 'Fp2-F7', 'Fp2-F8', 'Fp2-F3', 'Fp2-F4', ...
                'Fp2-T3', 'Fp2-T4', 'Fp2-P3', 'Fp2-P4', 'C3-C4', ...
                'C3-T5', 'C3-T6', 'C3-O1', 'C3-O2', 'C3-F7', ...
                'C3-F8', 'C3-F3', 'C3-F4', 'C3-T3', 'C3-T4', ...
                'C3-P3', 'C3-P4', 'C4-T5', 'C4-T6', 'C4-O1', ...
                'C4-O2', 'C4-F7', 'C4-F8', 'C4-F3', 'C4-F4', ...
                'C4-T3', 'C4-T4', 'C4-P3', 'C4-P4', 'T5-T6', ...
                'T5-O1', 'T5-O2', 'T5-F7', 'T5-F8', 'T5-F3', ...
                'T5-F4', 'T5-T3', 'T5-T4', 'T5-P3', 'T5-P4', ...
                'T6-O1', 'T6-O2', 'T6-F7', 'T6-F8', 'T6-F3', ...
                'T6-F4', 'T6-T3', 'T6-T4', 'T6-P3', 'T6-P4', ...
                'O1-O2', 'O1-F7', 'O1-F8', 'O1-F3', 'O1-F4', ...
                'O1-T3', 'O1-T4', 'O1-P3', 'O1-P4', 'O2-F7', ...
                'O2-F8', 'O2-F3', 'O2-F4', 'O2-T3', 'O2-T4', ...
                'O2-P3', 'O2-P4', 'F7-F8', 'F7-F3', 'F7-F4', ...
                'F7-T3', 'F7-T4', 'F7-P3', 'F7-P4', 'F8-F3', ...
                'F8-F4', 'F8-T3', 'F8-T4', 'F8-P3', 'F8-P4', ...
                'F3-F4', 'F3-T3', 'F3-T4', 'F3-P3', 'F3-P4', ...
                'F4-T3', 'F4-T4', 'F4-P3', 'F4-P4', 'T3-T4', ...
                'T3-P3', 'T3-P4', 'T4-P3', 'T4-P4', 'P3-P4'};
freqColors = [1, 0.75, 0.8; % pink
              0.5, 0, 0.5;  % purple
              0, 0, 1;      % blue
              0, 1, 0;      % green
              1, 1, 0;      % yellow
              1, 0.65, 0;   % orange
              1, 0, 0];     % red
freqNames = {'2Hz', '7Hz', '10Hz', '15Hz', '20Hz', '30Hz', '40Hz'};
numChannels = length(channels);
numChannelPairs = length(channelPairs);
numFrequencies = 7;
numTimePoints = 90000; % Always use 25 time points

% Electrode positions based on 10-20 system (in the same ***order*** as electrodeLabels)
electrodePos = [
    -0.125, 0.5;   % Fp1
    0.125, 0.5;    % Fp2
    -0.2, 0;       % C3
    0.2, 0;        % C4
    -0.3425, -0.27;% T5 (aka P7)
    0.3425, -0.27; % T6 (aka P8)
    -0.125, -0.5;  % O1
    0.125, -0.5;   % O2
    -0.3425, 0.27; % F7
    0.3425, 0.27;  % F8
    -0.18, 0.15;   % F3
    0.18, 0.15;    % F4
    -0.5, 0;       % T3
    0.5, 0;        % T4
    -0.18, -0.15;  % P3
    0.18, -0.15    % P4
];



% Define positions for the LED dots in a circle
theta = linspace(0, 2*pi, numFrequencies + 1);
theta(end) = []; % Remove last point to avoid duplication of the first point
ledPos = [cos(theta)', sin(theta)'];

% Initialize video writer
v = VideoWriter('phaseDiffVizSedated7to10HzRings.mp4', 'MPEG-4');
v.FrameRate = 25; % Adjust frame rate as needed
% Define expected size based on your observations
expectedSize = [1868, 3456];  % height, width
open(v);

% Fixed size for figure and axes
figurePosition = [0, 0, 3456, 1868]; % Define fixed figure size
axesPosition1 = [0.05, 0.1, 0.45, 0.8]; % Define fixed position for subplot 1
axesPosition2 = [0.55, 0.1, 0.4, 0.8]; % Define fixed position for subplot 2

% Loop through time points to update dynamic elements
for t = 1:numTimePoints
    % Create a new figure for each frame
    figh = figure('Position', figurePosition, 'Color', 'w'); % Set to specified frame size
    
    % Triangular correlation matrix
    ax1 = subplot('Position', axesPosition1); % Adjusted position to minimize white space
    hold(ax1, 'on');
    axis(ax1, 'equal');
    box(ax1, 'on'); % Keep the box around the plot
    set(ax1, 'XTick', 1:numChannels, 'XTickLabel', channels, 'YTick', 1:numChannels, 'YTickLabel', flip(channels), 'FontSize', 24);
    xtickangle(ax1, 45); % Set X-axis label rotation
    
    % Add text for box vs diamond
    annotation(figh, 'textbox', [0.38, 0.83, 0.1, 0.02], 'String', 'Connectivity:', 'FontSize', 16, 'EdgeColor', 'none');
    annotation(figh, 'textbox', [0.38, 0.81, 0.15, 0.02], 'String', 'Square = Left–right', 'FontSize', 16, 'EdgeColor', 'none');
    annotation(figh, 'textbox', [0.38, 0.79, 0.16, 0.02], 'String', 'Diamond = Front–back', 'FontSize', 16, 'EdgeColor', 'none');
    
    % Define the frequencies corresponding to each color
    frequencies = [2, 7, 10, 15, 20, 30, 40];

    % Add manual annotations for frequency-color mapping
    for i = 1:numFrequencies
        annotation(figh, 'textbox', [0.38, 0.77 - (i-1)*0.02, 0.15, 0.02], ...
                   'String', sprintf('%d Hz', frequencies(i)), ...
                   'FontSize', 16, 'EdgeColor', 'none', 'Color', freqColors(i, :));
    end
    
    % Title for the correlation matrix
    title(ax1, 'Triangular Correlation Matrix', 'FontSize', 24);

    % Define the noise gate threshold for triangular correlation matrix
    threshold = 30; % Phase alignment threshold in degrees

    % Plotting the triangular correlation matrix
    for i = 1:numChannels
        for j = 1:i-1
            chPairIdx = find(strcmp(channelPairs, sprintf('%s-%s', channels{i}, channels{j})));
            if isempty(chPairIdx)
                chPairIdx = find(strcmp(channelPairs, sprintf('%s-%s', channels{j}, channels{i})));
            end
            for freqIdx = 1:numFrequencies
                value = abs(data{freqIdx}(t, chPairIdx)); % Use absolute value to check alignment

                % Apply the noise gate threshold
                if value < threshold
                    % Calculate intensity based on alignment
                    intensity = (threshold - value) / threshold; % Scale intensity from 0 to 1
                    color = freqColors(freqIdx, :) * intensity; % Vary color intensity

                    % Position the dot in the triangular matrix
                    xPos = j + ledPos(freqIdx, 1) * 0.2;
                    yPos = (numChannels - i + 1) + ledPos(freqIdx, 2) * 0.2;

                    % Plot the dot in the triangular matrix
                    scatter(ax1, xPos, yPos, 100, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color, 'Marker', 'o');
                end
            end
            % Draw box or diamond around the LED circles
            if any(strcmp(leftRightComparisons(:,1), channels{i}) & strcmp(leftRightComparisons(:,2), channels{j}) | ...
                   strcmp(leftRightComparisons(:,2), channels{i}) & strcmp(leftRightComparisons(:,1), channels{j}))
                rectangle(ax1, 'Position', [j-0.4, numChannels-i+0.6, 0.8, 0.8], 'EdgeColor', 'k', 'LineWidth', 2);
            elseif any(strcmp(frontBackComparisons(:,1), channels{i}) & strcmp(frontBackComparisons(:,2), channels{j}) | ...
                       strcmp(frontBackComparisons(:,2), channels{i}) & strcmp(frontBackComparisons(:,1), channels{j}))
                plot(ax1, j, numChannels-i+1, 'kd', 'MarkerSize', 45, 'LineWidth', 3);
            end
        end
    end

    % Headplot with electrodes and connections
    ax2 = subplot('Position', axesPosition2); % Adjusted position and size
    hold(ax2, 'on');
    axis(ax2, 'equal');
    box(ax2, 'on'); % Keep the box around the plot
    set(ax2, 'XTick', [], 'YTick', []); % No tick marks, but keep axes

    % Plot head outline
    headX = [-0.75, 0.75, 0.75, -0.75, -0.75];
    headY = [1, 1, -1, -1, 1];
    plot(ax2, headX, headY, 'k');

    % Define the noise gate threshold for headplot
    threshold = 30; % Phase alignment threshold in degrees
    scaleFactor = 0.2;

    % Plot connections between the correct frequency dots
    for freqIdx = 1:numFrequencies
        for i = 1:numChannels
            for j = i+1:numChannels %j = 1:numChannels
                if i ~= j
                    % Get the index for the channel pair
                    chPairIdx = find(strcmp(channelPairs, sprintf('%s-%s', channels{i}, channels{j})));
                    if isempty(chPairIdx)
                        chPairIdx = find(strcmp(channelPairs, sprintf('%s-%s', channels{j}, channels{i})));
                    end
                    % Get the phase difference value
                    value = abs(data{freqIdx}(t, chPairIdx)); % Use absolute value to check alignment

                    % Apply the noise gate threshold
                    if value < threshold
                        % Calculate intensity based on alignment
                        intensity = (threshold - value) / threshold; % Scale intensity from 0 to 1
                        color = freqColors(freqIdx, :) * intensity; % Vary color intensity
                        
                        % Calculate the position on the ring for each electrode
                        ledX1 = electrodePos(i, 1) * 3 + ledPos(freqIdx, 1) * scaleFactor;
                        ledY1 = electrodePos(i, 2) * 3 + ledPos(freqIdx, 2) * scaleFactor;
                        ledX2 = electrodePos(j, 1) * 3 + ledPos(freqIdx, 1) * scaleFactor;
                        ledY2 = electrodePos(j, 2) * 3 + ledPos(freqIdx, 2) * scaleFactor;
                        
                        % Draw the line between the correct frequency positions
                        line(ax2, [ledX1, ledX2], [ledY1, ledY2], 'Color', color, 'LineWidth', 1 + intensity * 2);
                    end
                end
            end
        end
    end

    % Plot electrode labels on top
    for chIdx = 1:numChannels
        % Plot a white circle behind the label
        plot(ax2, electrodePos(chIdx, 1) * 3, electrodePos(chIdx, 2) * 3, 'wo', 'MarkerSize', 45, 'MarkerFaceColor', 'w');
        % Plot the electrode label in black
        text(ax2, electrodePos(chIdx, 1) * 3, electrodePos(chIdx, 2) * 3, channels{chIdx}, ...
            'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k');
    end

    hold(ax2, 'off');
    title(ax2, 'Headplot', 'FontSize', 24);

    % Capture frame for video
    frame = getframe(figh);

    % Check if the frame size matches the expected size
    frameSize = size(frame.cdata);  % Get the size of the frame [Height, Width, Channels]
    if all(frameSize(1:2) == expectedSize)  % Compare only Height and Width
        writeVideo(v, frame);
    else
        % Optionally, log a warning or handle this case
        warning('Frame size mismatch. Expected [%d, %d], but got [%d, %d]. Skipping this frame.', expectedSize(1), expectedSize(2), frameSize(1), frameSize(2));
    end

    % Close the figure to avoid memory issues
    close(figh);
end

% Close video writer
close(v);







%PhaseDensity!!
% Define color map with 120 unique colors
cmap = hsv(120);

% Initialize video writer
v = VideoWriter('phaseDensityPlotVideo.mp4', 'MPEG-4');
v.FrameRate = 1; % 1 frame per second for easy viewing, adjust as needed
open(v);

% Define the bins for phase lag (-180 to 180 degrees, with 5-degree steps)
bins = -180:5:180; % Bins from -180 to 180 degrees, with 5-degree steps
numBins = length(bins) - 1; % Number of bins

% Loop through each channel pair
for chPairIdx = 1:length(channelPairs)
    
    % Extract the relevant data for the current pair
    phaseData = phaseDiffDeg10Hz(:, chPairIdx);

    % Calculate the density of values in each bin
    density = histcounts(phaseData, bins);

    % Normalize the density for visualization
    density = density / max(density);

    % Create the plot
    theta = deg2rad(bins(1:end-1) + diff(bins)/2); % Convert bin centers to radians

    % Plot each segment as a line from the origin to a point on the circle
    figh = figure('Position', [0, 0, 800, 800], 'Color', 'w');
    hold on;
    
    % Plot the density segments
    for i = 1:numBins
        % Calculate the end point of the segment
        x = [0, density(i) * cos(theta(i))];
        y = [0, density(i) * sin(theta(i))];
        
        % Plot the segment using the color corresponding to the current channel pair
        plot(x, y, 'Color', cmap(chPairIdx, :), 'LineWidth', 2);
    end
    
    % Plot the circle boundary
    t = linspace(0, 2*pi, 100);
    plot(cos(t), sin(t), 'k--'); % Unit circle for reference

    % Formatting
    axis equal;
    title(sprintf('Phase Lag Density Plot for %s (10Hz)', channelPairs{chPairIdx}));
    xlabel('X');
    ylabel('Y');
    hold off;
    
    % Capture the frame for the video
    frame = getframe(figh);
    writeVideo(v, frame);
    
    % Close the figure to avoid memory issues
    close(figh);
end

% Close video writer
close(v);



%Phase Density Video!!
electrodePos = [
    -0.125, 0.5;   % Fp1
    0.125, 0.5;    % Fp2
    -0.2, 0;       % C3
    0.2, 0;        % C4
    -0.3425, -0.27;% T5 (aka P7)
    0.3425, -0.27; % T6 (aka P8)
    -0.125, -0.5;  % O1
    0.125, -0.5;   % O2
    -0.3425, 0.27; % F7
    0.3425, 0.27;  % F8
    -0.18, 0.15;   % F3
    0.18, 0.15;    % F4
    -0.5, 0;       % T3
    0.5, 0;        % T4
    -0.18, -0.15;  % P3
    0.18, -0.15    % P4
];

% Define color map for 7 frequencies
freqColors = [1, 0.75, 0.8; % pink
              0.5, 0, 0.5;  % purple
              0, 0, 1;      % blue
              0, 1, 0;      % green
              1, 1, 0;      % yellow
              1, 0.65, 0;   % orange
              1, 0, 0];     % red

% Frequency data
freqData = {phaseDiffDeg2Hz, phaseDiffDeg7Hz, phaseDiffDeg10Hz, ...
            phaseDiffDeg15Hz, phaseDiffDeg20Hz, phaseDiffDeg30Hz, ...
            phaseDiffDeg40Hz};

% Initialize video writer
v = VideoWriter('phaseDensityWithHeadplotVideo.mp4', 'MPEG-4');
v.FrameRate = 1; % 1 frame per second for easy viewing, adjust as needed
open(v);

% Define the bins for phase lag (-180 to 180 degrees, with 5-degree steps)
bins = -180:5:180; % Bins from -180 to 180 degrees, with 5-degree steps
numBins = length(bins) - 1; % Number of bins

% Define a larger angular offset for each frequency to avoid overlap
angularOffset = linspace(0, 3*pi/180, numFrequencies); % Larger offset in radians

% Loop through each channel pair
for chPairIdx = 1:length(channelPairs)
    
    % Create the figure
    figh = figure('Position', [0, 0, 1400, 800], 'Color', 'w');
    
    % Create a subplot for the phase density plot
    ax1 = subplot(1, 2, 1);
    hold(ax1, 'on');
    axis(ax1, 'equal');
    box(ax1, 'on');
    
    % Plot each frequency as a separate segment on the polar plot
    for freqIdx = 1:numFrequencies
        % Extract the relevant data for the current pair and frequency
        phaseData = freqData{freqIdx}(:, chPairIdx); 

        % Calculate the density of values in each bin
        density = histcounts(phaseData, bins);

        % Normalize the density for visualization
        density = density / max(density);

        % Apply the angular offset
        theta = deg2rad(bins(1:end-1) + diff(bins)/2) + angularOffset(freqIdx);

        % Plot the density segments
        for i = 1:numBins
            % Calculate the end point of the segment
            x = [0, density(i) * cos(theta(i))];
            y = [0, density(i) * sin(theta(i))];
            
            % Plot the segment using the color corresponding to the current frequency
            plot(ax1, x, y, 'Color', freqColors(freqIdx, :), 'LineWidth', 2);
        end
    end
    
    % Plot the circle boundary
    t = linspace(0, 2*pi, 100);
    plot(ax1, cos(t), sin(t), 'k--'); % Unit circle for reference
    
    % Add title
    title(ax1, sprintf('Phase Lag Density Plot for %s', channelPairs{chPairIdx}));
    
    % Create the headplot subplot
    ax2 = subplot(1, 2, 2);
    hold(ax2, 'on');
    axis(ax2, 'equal');
    box(ax2, 'on');
    set(ax2, 'XTick', [], 'YTick', []); % No tick marks, but keep axes
    
    % Plot head outline
    headX = [-0.75, 0.75, 0.75, -0.75, -0.75];
    headY = [1, 1, -1, -1, 1];
    plot(ax2, headX, headY, 'k');
    
    % Extract the channel labels for the current pair
    chPair = strsplit(channelPairs{chPairIdx}, '-');
    ch1 = find(strcmp(channels, chPair{1}));
    ch2 = find(strcmp(channels, chPair{2}));

    % Plot the line connecting the two electrodes
    ledX1 = electrodePos(ch1, 1) * 3;
    ledY1 = electrodePos(ch1, 2) * 3;
    ledX2 = electrodePos(ch2, 1) * 3;
    ledY2 = electrodePos(ch2, 2) * 3;
    plot(ax2, [ledX1, ledX2], [ledY1, ledY2], 'k-', 'LineWidth', 2);

    % Plot the electrode labels
    for chIdx = 1:numChannels
        % Plot a white circle behind the label
        plot(ax2, electrodePos(chIdx, 1) * 3, electrodePos(chIdx, 2) * 3, 'wo', 'MarkerSize', 45, 'MarkerFaceColor', 'w');
        % Plot the electrode label in black
        text(ax2, electrodePos(chIdx, 1) * 3, electrodePos(chIdx, 2) * 3, channels{chIdx}, ...
            'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k');
    end
    
    title(ax2, sprintf('Headplot for %s', channelPairs{chPairIdx}));
    
    % Manually create the legend with correct colors
    legendLabels = gobjects(1, numFrequencies);
    for freqIdx = 1:numFrequencies
        legendLabels(freqIdx) = plot(nan, nan, 'Color', freqColors(freqIdx, :), 'LineWidth', 2);
    end
    legend(ax1, legendLabels, freqNames, 'Location', 'northeastoutside');
    
    % Capture the frame for the video
    frame = getframe(figh);
    writeVideo(v, frame);

    % Close the figure to avoid memory issues
    close(figh);
end

% Close video writer
close(v);






%Mean phase lag value at any given frequency at Fp1:F3.
% Index for Fp1-F3 (which is item 10 in the channelPairs list)
chPairIdx = 10;

% Frequency data
freqData = {phaseDiffDeg2Hz, phaseDiffDeg7Hz, phaseDiffDeg10Hz, ...
            phaseDiffDeg15Hz, phaseDiffDeg20Hz, phaseDiffDeg30Hz, ...
            phaseDiffDeg40Hz};

% Initialize array to store the average phase lag for each frequency
averagePhaseLag = zeros(1, numFrequencies);

% Loop through each frequency and compute the average phase lag for Fp1-F3
for freqIdx = 1:numFrequencies
    % Extract the phase lag data for Fp1-F3 (in degrees)
    phaseData = freqData{freqIdx}(:, chPairIdx);
    
    % Convert phase data from degrees to radians
    phaseDataRad = deg2rad(phaseData);
    
    % Compute the circular mean of the phase data
    avgPhaseRad = angle(mean(exp(1i * phaseDataRad)));
    
    % Convert the result back to degrees
    averagePhaseLag(freqIdx) = rad2deg(avgPhaseRad);
end

% Display the average phase lag for each frequency
for freqIdx = 1:numFrequencies
    fprintf('Average phase lag for %s at %d Hz: %.2f degrees\n', channelPairs{chPairIdx}, frequencies(freqIdx), averagePhaseLag(freqIdx));
end


%Mean phase lag and amplitude coalescence at each freq. Fp1-F3.
% Index for Fp1-F3 (item 10 in the channelPairs list)
chPairIdx = 10;

% Frequency data
freqData = {phaseDiffDeg2Hz, phaseDiffDeg7Hz, phaseDiffDeg10Hz, ...
            phaseDiffDeg15Hz, phaseDiffDeg20Hz, phaseDiffDeg30Hz, ...
            phaseDiffDeg40Hz};

% Define frequencies and their colors
frequencies = [2, 7, 10, 15, 20, 30, 40];
freqColors = [1, 0.75, 0.8; % pink
              0.5, 0, 0.5;  % purple
              0, 0, 1;      % blue
              0, 1, 0;      % green
              1, 1, 0;      % yellow
              1, 0.65, 0;   % orange
              1, 0, 0];     % red

% Initialize arrays for phase directions and amplitudes
averagePhaseLag = zeros(1, numFrequencies);  % Average phase lag (degrees)
vectorStrength = zeros(1, numFrequencies);   % Vector strength (amplitude)

% Loop through each frequency and compute the average phase lag and vector strength
for freqIdx = 1:numFrequencies
    % Extract the phase lag data for Fp1-F3 (in degrees)
    phaseData = freqData{freqIdx}(:, chPairIdx);
    
    % Convert phase data from degrees to radians
    phaseDataRad = deg2rad(phaseData);
    
    % Compute the circular mean of the phase data
    avgPhaseRad = angle(mean(exp(1i * phaseDataRad)));
    
    % Compute vector strength (coalescence amplitude)
    vectorStrength(freqIdx) = abs(mean(exp(1i * phaseDataRad)));
    
    % Convert the result back to degrees
    averagePhaseLag(freqIdx) = rad2deg(avgPhaseRad);
end

% Plot the average phase lag as a polar plot
figure;
polarAxes = polaraxes;
hold on;

% Loop to plot each frequency
for freqIdx = 1:numFrequencies
    % Convert the average phase lag to radians for polar plot
    avgPhaseRad = deg2rad(averagePhaseLag(freqIdx));
    
    % Plot a vector for each frequency
    polarplot([0 avgPhaseRad], [0 vectorStrength(freqIdx)], 'Color', freqColors(freqIdx, :), 'LineWidth', 2);
end

% Customize the polar plot
title('Average Phase Lag for Fp1-F3 Across Frequencies');
set(polarAxes, 'ThetaZeroLocation', 'right'); % Set 0 degrees to the right
set(polarAxes, 'ThetaDir', 'counterclockwise'); % Set angle increase direction
legend({'2 Hz', '7 Hz', '10 Hz', '15 Hz', '20 Hz', '30 Hz', '40 Hz'}, 'Location', 'southoutside');

hold off;







% Render a video of bandpass time series and red envelope. 2160x2160 pixels.

% Data details
samplingRate = 125; % 125 Hz sampling rate
windowSize = 625; % 5 seconds window (625 time points)
stepSize = 5; % Advance 5 time points per frame
frameCount = 89998; % Number of frames to display (for now)

% Load your data (assuming the CSV is stored as complex numbers)
complexNumbers2Hz = csvread('complexNumbers2Hz.csv');
complexNumbers2Hz = reshape(complexNumbers2Hz, [], 16); 
complexNumbers5Hz = csvread('complexNumbers5Hz.csv');
complexNumbers5Hz = reshape(complexNumbers5Hz, [], 16); 
complexNumbers7Hz = csvread('complexNumbers7Hz.csv');
complexNumbers7Hz = reshape(complexNumbers7Hz, [], 16); 

% Correct order of channels from top to bottom
electrodeLabels = {'Fp1', 'Fp2', 'F7', 'F3', 'F4', 'F8', 'T3', 'C3', 'C4', 'T4', ...
                   'T5', 'P3', 'P4', 'T6', 'O1', 'O2'};
electrodeLabelsBackwards = {'O2', 'O1', 'T6', 'P4', 'P3', 'T5', 'T4', 'C4', 'C3', 'T3', 'F8', 'F4', 'F3', 'F7', 'Fp2', 'Fp1'};

% Reorder the data according to the new electrode label order
reorderedData2 = complexNumbers2Hz(:, [8, 7, 6, 16, 15, 5, 14, 4, 3, 13, ...
                                       10, 12, 11, 9, 2, 1]);
reorderedData5 = complexNumbers2Hz(:, [8, 7, 6, 16, 15, 5, 14, 4, 3, 13, ...
                                       10, 12, 11, 9, 2, 1]);
reorderedData7 = complexNumbers2Hz(:, [8, 7, 6, 16, 15, 5, 14, 4, 3, 13, ...
                                       10, 12, 11, 9, 2, 1]);
% Initialize figure for the video (starting square at 1000x1000)
figh = figure('Position', [100, 100, 1000, 1000], 'Color', 'w');
axis equal;
axis off; % Hide axes

% Setup video writer
videoFileName = fullfile(getenv('HOME'), 'Desktop', 'EEG_Bandpass2Hz_Envelope.mp4'); % Save video to Desktop
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 25; % Set the frame rate
open(v);

% Initialize time point and empty plot (cursor starts at the right)
currentDataIndex = 1; % Start at the first data point

% Setup time window loop
for frameNum = 1:frameCount
    clf; % Clear figure

    % Left column: Time series + envelope
    hold on;
    for ch = 1:16
        % Extract current window of data for the channel
        if currentDataIndex >= windowSize
            currentData = reorderedData(currentDataIndex-windowSize+1:currentDataIndex, ch);
        else
            currentData = [zeros(windowSize-currentDataIndex, 1); reorderedData(1:currentDataIndex, ch)];
        end

        % Shift the channel vertically
        verticalOffset = 30 * (ch - 1.5); % Adjust spacing between channels

        % Plot the real part of the signal (flowing leftwards)
        plot(linspace(0, 5, windowSize), real(currentData) + verticalOffset, 'b', 'LineWidth', 2);

        % Plot the envelope (absolute value of complex numbers)
        plot(linspace(0, 5, windowSize), abs(currentData) + verticalOffset, 'r', 'LineWidth', 2);

        % Add channel labels on the left (only once)
        
            text(-0.5, verticalOffset, electrodeLabelsBackwards{ch}, 'FontSize', 14, 'Color', 'k');
        
    end

    % Calculate current time (in seconds)
    currentTime = (currentDataIndex - 1) / samplingRate;

    % Convert to hours, minutes, and seconds format
    hours = floor(currentTime / 3600);
    minutes = floor(mod(currentTime, 3600) / 60);
    seconds = mod(currentTime, 60);

    % Add recording time below the plot
    timeString = sprintf('Recording time: %02d:%02d:%04.1f s', hours, minutes, seconds);
    text(2.5, -70, timeString, 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', 'k');

    % Remove y-axis label
    set(gca, 'ytick', []);
    
    title('EEG Time Series (2 Hz Bandpass) with Envelope');
    axis([0 5 -50 450]); % Keep axis limits constant, allowing room for labels
    set(gca, 'FontSize', 14); % Set font size for better video readability
    drawnow;

    % Capture frame and write to video
    frame = getframe(figh);
    writeVideo(v, frame);
    
    % Move the cursor to the right and let data scroll leftwards
    currentDataIndex = currentDataIndex + stepSize;

    % Check if we reach the end of the data
    if currentDataIndex > size(reorderedData, 1)
        break;
    end
end

% Close video writer
close(v);
close(figh); % Close the figure

disp('Video rendering complete. Check your Desktop.');






%3 at once!
% Data details
samplingRate = 125; % 125 Hz sampling rate
windowSize = 625; % 5 seconds window (625 time points)
stepSize = 5; % Advance 5 time points per frame
frameCount = 89998; % Number of frames to display (for now)

% Load your data (assuming the CSV is stored as complex numbers)
complexNumbers2Hz = csvread('complexNumbers2Hz.csv');
complexNumbers2Hz = reshape(complexNumbers2Hz, [], 16); 
complexNumbers5Hz = csvread('complexNumbers5Hz.csv');
complexNumbers5Hz = reshape(complexNumbers5Hz, [], 16); 
complexNumbers7Hz = csvread('complexNumbers7Hz.csv');
complexNumbers7Hz = reshape(complexNumbers7Hz, [], 16); 

% Correct order of channels from top to bottom
electrodeLabels = {'Fp1', 'Fp2', 'F7', 'F3', 'F4', 'F8', 'T3', 'C3', 'C4', 'T4', ...
                   'T5', 'P3', 'P4', 'T6', 'O1', 'O2'};
electrodeLabelsBackwards = {'O2', 'O1', 'T6', 'P4', 'P3', 'T5', 'T4', 'C4', 'C3', 'T3', 'F8', 'F4', 'F3', 'F7', 'Fp2', 'Fp1'};

% Reorder the data according to the new electrode label order
reorderedData2 = complexNumbers2Hz(:, [8, 7, 6, 16, 15, 5, 14, 4, 3, 13, ...
                                       10, 12, 11, 9, 2, 1]);
reorderedData5 = complexNumbers5Hz(:, [8, 7, 6, 16, 15, 5, 14, 4, 3, 13, ...
                                       10, 12, 11, 9, 2, 1]);
reorderedData7 = complexNumbers7Hz(:, [8, 7, 6, 16, 15, 5, 14, 4, 3, 13, ...
                                       10, 12, 11, 9, 2, 1]);

% Initialize figure for the video (starting square at 1000x1000)
figh = figure('Position', [100, 100, 1000, 1000], 'Color', 'w');
axis equal;
axis off; % Hide axes

% Setup video writer
videoFileName = fullfile(getenv('HOME'), 'Desktop', 'EEG_Bandpass2Hz_5Hz_7Hz_Envelope.mp4'); % Save video to Desktop
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 25; % Set the frame rate
open(v);

% Initialize time point and empty plot (cursor starts at the right)
currentDataIndex = 1; % Start at the first data point

% Setup time window loop
for frameNum = 1:frameCount
    clf; % Clear figure

    % Left column: Time series + envelope
    hold on;
    for ch = 1:16
        % Extract current window of data for each bandpass channel
        if currentDataIndex >= windowSize
            currentData2 = reorderedData2(currentDataIndex-windowSize+1:currentDataIndex, ch);
            currentData5 = reorderedData5(currentDataIndex-windowSize+1:currentDataIndex, ch);
            currentData7 = reorderedData7(currentDataIndex-windowSize+1:currentDataIndex, ch);
        else
            currentData2 = [zeros(windowSize-currentDataIndex, 1); reorderedData2(1:currentDataIndex, ch)];
            currentData5 = [zeros(windowSize-currentDataIndex, 1); reorderedData5(1:currentDataIndex, ch)];
            currentData7 = [zeros(windowSize-currentDataIndex, 1); reorderedData7(1:currentDataIndex, ch)];
        end

        % Shift the channel vertically
        verticalOffset = 30 * (ch - 1.5); % Adjust spacing between channels

scalingFactor = 1.0;  % Scale down the amplitude to make the plot look better

% Plot the 2 Hz bandpass signal and envelope
plot(linspace(0, 5, windowSize), real(currentData2) * scalingFactor + verticalOffset, 'Color', [0 0 0.5], 'LineWidth', 3.0); % Dark blue for 2 Hz
plot(linspace(0, 5, windowSize), abs(currentData2) * scalingFactor + verticalOffset, 'Color', [0.68 0.85 0.9], 'LineWidth', 2.0); % Light blue for 2 Hz envelope

% Plot the 5 Hz bandpass signal and envelope
plot(linspace(0, 5, windowSize), real(currentData5) * scalingFactor + verticalOffset, 'Color', [0 0.5 0], 'LineWidth', 3.0); % Dark green for 5 Hz
plot(linspace(0, 5, windowSize), abs(currentData5) * scalingFactor + verticalOffset, 'Color', [0.56 0.93 0.56], 'LineWidth', 2.0); % Light green for 5 Hz envelope

% Plot the 7 Hz bandpass signal and envelope
plot(linspace(0, 5, windowSize), real(currentData7) * scalingFactor + verticalOffset, 'Color', [0.5 0 0], 'LineWidth', 3.0); % Dark red for 7 Hz
plot(linspace(0, 5, windowSize), abs(currentData7) * scalingFactor + verticalOffset, 'Color', [1 0.71 0.76], 'LineWidth', 2.0); % Pink for 7 Hz envelope

        % Add channel labels on the left (only once)
        if frameNum == 1
            text(-0.5, verticalOffset, electrodeLabelsBackwards{ch}, 'FontSize', 14, 'Color', 'k');
        end
    end

    % Calculate current time (in seconds)
    currentTime = (currentDataIndex - 1) / samplingRate;

    % Convert to hours, minutes, and seconds format
    hours = floor(currentTime / 3600);
    minutes = floor(mod(currentTime, 3600) / 60);
    seconds = mod(currentTime, 60);

    % Add recording time below the plot
    timeString = sprintf('Recording time: %02d:%02d:%04.1f s', hours, minutes, seconds);
    text(2.5, -70, timeString, 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', 'k');

    % Remove y-axis label
    set(gca, 'ytick', []);

    % Add title indicating bandpass colors
    title({'EEG Time Series for 2 Hz (Blue/LightBlueEnv), 5 Hz (Green/LightGreenEnv), 7 Hz (Red/PinkEnv)'});
    axis([0 5 -50 450]); % Keep axis limits constant, allowing room for labels
    set(gca, 'FontSize', 14); % Set font size for better video readability
    drawnow;

    % Capture frame and write to video
    frame = getframe(figh);
    writeVideo(v, frame);
    
    % Move the cursor to the right and let data scroll leftwards
    currentDataIndex = currentDataIndex + stepSize;

    % Check if we reach the end of the data
    if currentDataIndex > size(reorderedData2, 1)
        break;
    end
end

% Close video writer
close(v);
close(figh); % Close the figure

disp('Video rendering complete. Check your Desktop.');





% Now, a video of helices with fading vectors, dots, and present-time vector for 2 Hz, 5 Hz, 7 Hz, and 10 Hz!

% Data details
samplingRate = 125; % 125 Hz sampling rate
stepSize = 5; % Time point step size
frameCount = 15000; % Number of frames to display (adjust as needed)
fadeLength = 100; % Length of time fade (how many points to show fading)
windowSize = 625; % 5 seconds window

% Load your data (assuming the CSV is stored as complex numbers)
complexNumbers2Hz = csvread('complexNumbers2Hz.csv');
complexNumbers2Hz = reshape(complexNumbers2Hz, [], 16); 
complexNumbers5Hz = csvread('complexNumbers5Hz.csv');
complexNumbers5Hz = reshape(complexNumbers5Hz, [], 16);
complexNumbers7Hz = csvread('complexNumbers7Hz.csv');
complexNumbers7Hz = reshape(complexNumbers7Hz, [], 16);
complexNumbers10Hz = csvread('complexNumbers10Hz.csv');
complexNumbers10Hz = reshape(complexNumbers10Hz, [], 16);

% Use the 6th column for T6 and the 2nd column for Fp2
T6_data_2Hz = complexNumbers2Hz(:, 6);
Fp2_data_2Hz = complexNumbers2Hz(:, 2);
T6_data_5Hz = complexNumbers5Hz(:, 6);
Fp2_data_5Hz = complexNumbers5Hz(:, 2);
T6_data_7Hz = complexNumbers7Hz(:, 6);
Fp2_data_7Hz = complexNumbers7Hz(:, 2);
T6_data_10Hz = complexNumbers10Hz(:, 6);
Fp2_data_10Hz = complexNumbers10Hz(:, 2);

% Initialize figure for the video (doubled in width)
figh = figure('Position', [100, 100, 2000, 1000], 'Color', 'w');

% Setup video writer
videoFileName = fullfile(getenv('HOME'), 'Desktop', 'T6_and_Fp2_ComplexPlane_MultipleHz_withFadingVectors.mp4'); % Save video to Desktop
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 25; % Set the frame rate
open(v);

% Initialize time point and plot
currentDataIndex = 1; % Start at the first data point
z_offset = 0; % Offset for time fade in z-axis

% Setup time window loop
for frameNum = 1:frameCount
    clf; % Clear figure

    % Create subplot for T6
    subplot(1, 2, 1);
    hold on;

        % Calculate current time (in seconds)
    currentTime = frameNum / samplingRate;

    % Convert to hours, minutes, and seconds format
    hours = floor(currentTime / 3600);
    minutes = floor(mod(currentTime, 3600) / 60);
    seconds = mod(currentTime, 60);

    % Add recording time below the plot
    timeString = sprintf('Recording time: %02d:%02d:%04.1f s', hours, minutes, seconds);
    text(2.5, -70, timeString, 'FontSize', 14, 'HorizontalAlignment', 'center', 'Color', 'k');


    % Plot black reference line at real=0, imag=0 for T6
    plot3([0 0], [0 0], [0 fadeLength], 'k-', 'LineWidth', 2);

    % Plot 2 Hz for T6
    if currentDataIndex >= fadeLength
        currentWindow2Hz = T6_data_2Hz(currentDataIndex-fadeLength+1:currentDataIndex);
    else
        currentWindow2Hz = [zeros(fadeLength-currentDataIndex, 1); T6_data_2Hz(1:currentDataIndex)];
    end
    realParts2Hz = real(currentWindow2Hz);
    imagParts2Hz = imag(currentWindow2Hz);
    z_vals2Hz = z_offset + (1:length(currentWindow2Hz)); % Reverse the time axis
    for i = 1:length(currentWindow2Hz)
        opacity = (i / fadeLength);
        fadeColor = [0, 0, 1] * opacity + [1, 1, 1] * (1 - opacity); % Blue fade
        plot3([0 realParts2Hz(i)], [0 imagParts2Hz(i)], [z_vals2Hz(i) z_vals2Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
        if i > 1
            plot3([realParts2Hz(i-1) realParts2Hz(i)], [imagParts2Hz(i-1) imagParts2Hz(i)], [z_vals2Hz(i-1) z_vals2Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
        end
        plot3(realParts2Hz(i), imagParts2Hz(i), z_vals2Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
    end

    % Repeat the above block for other frequencies (5 Hz, 7 Hz, 10 Hz) for T6
    % Plot the 5 Hz frequency in green, 7 Hz in gold, 10 Hz in red (same as original code)
    % ...

    % Plot 5 Hz for T6
if currentDataIndex >= fadeLength
    currentWindow5Hz = T6_data_5Hz(currentDataIndex-fadeLength+1:currentDataIndex);
else
    currentWindow5Hz = [zeros(fadeLength-currentDataIndex, 1); T6_data_5Hz(1:currentDataIndex)];
end
realParts5Hz = real(currentWindow5Hz);
imagParts5Hz = imag(currentWindow5Hz);
z_vals5Hz = z_offset + (1:length(currentWindow5Hz)); % Reverse the time axis
for i = 1:length(currentWindow5Hz)
    opacity = (i / fadeLength);
    fadeColor = [0, 1, 0] * opacity + [1, 1, 1] * (1 - opacity); % Green fade for 5 Hz
    plot3([0 realParts5Hz(i)], [0 imagParts5Hz(i)], [z_vals5Hz(i) z_vals5Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    if i > 1
        plot3([realParts5Hz(i-1) realParts5Hz(i)], [imagParts5Hz(i-1) imagParts5Hz(i)], [z_vals5Hz(i-1) z_vals5Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    end
    plot3(realParts5Hz(i), imagParts5Hz(i), z_vals5Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
end


% Plot 7 Hz for T6
if currentDataIndex >= fadeLength
    currentWindow7Hz = T6_data_7Hz(currentDataIndex-fadeLength+1:currentDataIndex);
else
    currentWindow7Hz = [zeros(fadeLength-currentDataIndex, 1); T6_data_7Hz(1:currentDataIndex)];
end
realParts7Hz = real(currentWindow7Hz);
imagParts7Hz = imag(currentWindow7Hz);
z_vals7Hz = z_offset + (1:length(currentWindow7Hz)); % Reverse the time axis
for i = 1:length(currentWindow7Hz)
    opacity = (i / fadeLength);
    fadeColor = [1, 0.84, 0] * opacity + [1, 1, 1] * (1 - opacity); % Gold fade for 7 Hz
    plot3([0 realParts7Hz(i)], [0 imagParts7Hz(i)], [z_vals7Hz(i) z_vals7Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    if i > 1
        plot3([realParts7Hz(i-1) realParts7Hz(i)], [imagParts7Hz(i-1) imagParts7Hz(i)], [z_vals7Hz(i-1) z_vals7Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    end
    plot3(realParts7Hz(i), imagParts7Hz(i), z_vals7Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
end

% Plot 10 Hz for T6
if currentDataIndex >= fadeLength
    currentWindow10Hz = T6_data_10Hz(currentDataIndex-fadeLength+1:currentDataIndex);
else
    currentWindow10Hz = [zeros(fadeLength-currentDataIndex, 1); T6_data_10Hz(1:currentDataIndex)];
end
realParts10Hz = real(currentWindow10Hz);
imagParts10Hz = imag(currentWindow10Hz);
z_vals10Hz = z_offset + (1:length(currentWindow10Hz)); % Reverse the time axis
for i = 1:length(currentWindow10Hz)
    opacity = (i / fadeLength);
    fadeColor = [1, 0, 0] * opacity + [1, 1, 1] * (1 - opacity); % Red fade for 10 Hz
    plot3([0 realParts10Hz(i)], [0 imagParts10Hz(i)], [z_vals10Hz(i) z_vals10Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    if i > 1
        plot3([realParts10Hz(i-1) realParts10Hz(i)], [imagParts10Hz(i-1) imagParts10Hz(i)], [z_vals10Hz(i-1) z_vals10Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    end
    plot3(realParts10Hz(i), imagParts10Hz(i), z_vals10Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
end


    xlabel('Imaginary');
    ylabel('Real');
    zlabel('Time fade');
    title('T6 Complex Plane for 2 Hz, 5 Hz, 7 Hz, 10 Hz with Fading Vectors');
    grid on;
    view([0.1 0.1 0.2]); % Adjust this to your desired viewing angle
    axis([-20 20 -20 20 0 fadeLength]); % Reverse the time axis (time increasing upwards)

    % Create subplot for Fp2
    subplot(1, 2, 2);
    hold on;

    % Plot black reference line at real=0, imag=0 for Fp2
    plot3([0 0], [0 0], [0 fadeLength], 'k-', 'LineWidth', 2);

    % Plot 2 Hz for Fp2
    if currentDataIndex >= fadeLength
        currentWindow2Hz = Fp2_data_2Hz(currentDataIndex-fadeLength+1:currentDataIndex);
    else
        currentWindow2Hz = [zeros(fadeLength-currentDataIndex, 1); Fp2_data_2Hz(1:currentDataIndex)];
    end
    realParts2Hz = real(currentWindow2Hz);
    imagParts2Hz = imag(currentWindow2Hz);
    z_vals2Hz = z_offset + (1:length(currentWindow2Hz)); % Reverse the time axis
    for i = 1:length(currentWindow2Hz)
        opacity = (i / fadeLength);
        fadeColor = [0, 0, 1] * opacity + [1, 1, 1] * (1 - opacity); % Blue fade
        plot3([0 realParts2Hz(i)], [0 imagParts2Hz(i)], [z_vals2Hz(i) z_vals2Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
        if i > 1
            plot3([realParts2Hz(i-1) realParts2Hz(i)], [imagParts2Hz(i-1) imagParts2Hz(i)], [z_vals2Hz(i-1) z_vals2Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
        end
        plot3(realParts2Hz(i), imagParts2Hz(i), z_vals2Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
    end

    % Repeat the above block for other frequencies (5 Hz, 7 Hz, 10 Hz) for Fp2
    % Plot the 5 Hz frequency in green, 7 Hz in gold, 10 Hz in red (same as original code)
    % ...

% Plot 5 Hz for Fp2
if currentDataIndex >= fadeLength
    currentWindow5Hz = Fp2_data_5Hz(currentDataIndex-fadeLength+1:currentDataIndex);
else
    currentWindow5Hz = [zeros(fadeLength-currentDataIndex, 1); Fp2_data_5Hz(1:currentDataIndex)];
end
realParts5Hz = real(currentWindow5Hz);
imagParts5Hz = imag(currentWindow5Hz);
z_vals5Hz = z_offset + (1:length(currentWindow5Hz)); % Reverse the time axis
for i = 1:length(currentWindow5Hz)
    opacity = (i / fadeLength);
    fadeColor = [0, 1, 0] * opacity + [1, 1, 1] * (1 - opacity); % Green fade for 5 Hz
    plot3([0 realParts5Hz(i)], [0 imagParts5Hz(i)], [z_vals5Hz(i) z_vals5Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    if i > 1
        plot3([realParts5Hz(i-1) realParts5Hz(i)], [imagParts5Hz(i-1) imagParts5Hz(i)], [z_vals5Hz(i-1) z_vals5Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    end
    plot3(realParts5Hz(i), imagParts5Hz(i), z_vals5Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
end

% Plot 7 Hz for Fp2
if currentDataIndex >= fadeLength
    currentWindow7Hz = Fp2_data_7Hz(currentDataIndex-fadeLength+1:currentDataIndex);
else
    currentWindow7Hz = [zeros(fadeLength-currentDataIndex, 1); Fp2_data_7Hz(1:currentDataIndex)];
end
realParts7Hz = real(currentWindow7Hz);
imagParts7Hz = imag(currentWindow7Hz);
z_vals7Hz = z_offset + (1:length(currentWindow7Hz)); % Reverse the time axis
for i = 1:length(currentWindow7Hz)
    opacity = (i / fadeLength);
    fadeColor = [1, 0.84, 0] * opacity + [1, 1, 1] * (1 - opacity); % Gold fade for 7 Hz
    plot3([0 realParts7Hz(i)], [0 imagParts7Hz(i)], [z_vals7Hz(i) z_vals7Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    if i > 1
        plot3([realParts7Hz(i-1) realParts7Hz(i)], [imagParts7Hz(i-1) imagParts7Hz(i)], [z_vals7Hz(i-1) z_vals7Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    end
    plot3(realParts7Hz(i), imagParts7Hz(i), z_vals7Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
end

% Plot 10 Hz for Fp2
if currentDataIndex >= fadeLength
    currentWindow10Hz = Fp2_data_10Hz(currentDataIndex-fadeLength+1:currentDataIndex);
else
    currentWindow10Hz = [zeros(fadeLength-currentDataIndex, 1); Fp2_data_10Hz(1:currentDataIndex)];
end
realParts10Hz = real(currentWindow10Hz);
imagParts10Hz = imag(currentWindow10Hz);
z_vals10Hz = z_offset + (1:length(currentWindow10Hz)); % Reverse the time axis
for i = 1:length(currentWindow10Hz)
    opacity = (i / fadeLength);
    fadeColor = [1, 0, 0] * opacity + [1, 1, 1] * (1 - opacity); % Red fade for 10 Hz
    plot3([0 realParts10Hz(i)], [0 imagParts10Hz(i)], [z_vals10Hz(i) z_vals10Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    if i > 1
        plot3([realParts10Hz(i-1) realParts10Hz(i)], [imagParts10Hz(i-1) imagParts10Hz(i)], [z_vals10Hz(i-1) z_vals10Hz(i)], 'Color', fadeColor, 'LineWidth', 1);
    end
    plot3(realParts10Hz(i), imagParts10Hz(i), z_vals10Hz(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
end


    xlabel('Imaginary');
    ylabel('Real');
    zlabel('Time fade');
    title('Fp2 Complex Plane for 2 Hz, 5 Hz, 7 Hz, 10 Hz with Fading Vectors');
    grid on;
    view([0.1 0.1 0.2]); % Adjust this to your desired viewing angle
    axis([-20 20 -20 20 0 fadeLength]); % Reverse the time axis (time increasing upwards)

    % Update the current data index
    currentDataIndex = currentDataIndex + stepSize;

    % Capture the frame
    frame = getframe(figh);
    writeVideo(v, frame);

    % Check if we reach the end of the data
    if currentDataIndex > size(T6_data_2Hz, 1)
        break;
    end
end

% Close video writer
close(v);
close(figh); % Close the figure






%Complex as in MultiFrek and also complex as in real+imag looking at
%complex EEG!!!

% Number of time points and electrodes
numTimePoints = 449999;
numChannels = 16;

% Load your complex EEG data for each frequency band
complexNumbers2Hz = csvread('complexNumbers2Hz.csv');
complexNumbers2Hz = reshape(complexNumbers2Hz, [], numChannels); 

complexNumbers5Hz = csvread('complexNumbers5Hz.csv');
complexNumbers5Hz = reshape(complexNumbers5Hz, [], numChannels); 

complexNumbers7Hz = csvread('complexNumbers7Hz.csv');
complexNumbers7Hz = reshape(complexNumbers7Hz, [], numChannels); 

complexNumbers10Hz = csvread('complexNumbers10Hz.csv');
complexNumbers10Hz = reshape(complexNumbers10Hz, [], numChannels); 

complexNumbers15Hz = csvread('complexNumbers15Hz.csv');
complexNumbers15Hz = reshape(complexNumbers15Hz, [], numChannels); 

complexNumbers20Hz = csvread('complexNumbers20Hz.csv');
complexNumbers20Hz = reshape(complexNumbers20Hz, [], numChannels); 

complexNumbers30Hz = csvread('complexNumbers30Hz.csv');
complexNumbers30Hz = reshape(complexNumbers30Hz, [], numChannels); 

complexNumbers40Hz = csvread('complexNumbers40Hz.csv');
complexNumbers40Hz = reshape(complexNumbers40Hz, [], numChannels); 

% Initialize matrix for the combined complex values (multiFrekComplex)
multiFrekComplex = zeros(numTimePoints, numChannels);

% Sum the real and imaginary parts across all frequencies for each time point and channel
for t = 1:numTimePoints
    for ch = 1:numChannels
        % Sum the real parts of each frequency at this time point
        realSum = real(complexNumbers2Hz(t, ch)) + real(complexNumbers5Hz(t, ch)) + ...
                  real(complexNumbers7Hz(t, ch)) + real(complexNumbers10Hz(t, ch)) + ...
                  real(complexNumbers15Hz(t, ch)) + real(complexNumbers20Hz(t, ch)) + ...
                  real(complexNumbers30Hz(t, ch)) + real(complexNumbers40Hz(t, ch));
              
        % Sum the imaginary parts of each frequency at this time point
        imagSum = imag(complexNumbers2Hz(t, ch)) + imag(complexNumbers5Hz(t, ch)) + ...
                  imag(complexNumbers7Hz(t, ch)) + imag(complexNumbers10Hz(t, ch)) + ...
                  imag(complexNumbers15Hz(t, ch)) + imag(complexNumbers20Hz(t, ch)) + ...
                  imag(complexNumbers30Hz(t, ch)) + imag(complexNumbers40Hz(t, ch));
              
        % Store the summed complex number for this time point and channel
        multiFrekComplex(t, ch) = realSum + 1i * imagSum;
    end
end

% Now multiFrekComplex is a 449999x16 matrix, containing the summed complex values
disp('multiFrekComplex computation complete.');

% You can save this result for later use
csvwrite('multiFrekComplex.csv', multiFrekComplex);







% Now, a video of helices with fading vectors, dots, and present-time vector for multiFrekComplex!

% Data details
samplingRate = 125; % 125 Hz sampling rate
stepSize = 5; % Time point step size
frameCount = 500; % Number of frames to display (adjust as needed)
fadeLength = 100; % Length of time fade (how many points to show fading)
windowSize = 625; % 5 seconds window

% Assuming multiFrekComplex is already calculated and is 449999 time points x 16 channels
% For T6 (column 6) and Fp2 (column 2), we extract the complex data:
T6_data = multiFrekComplex(:, 6); 
Fp2_data = multiFrekComplex(:, 2);

% Initialize figure for the video (doubled in width)
figh = figure('Position', [100, 100, 2000, 1000], 'Color', 'w');

% Setup video writer
videoFileName = fullfile(getenv('HOME'), 'Desktop', 'T6_and_Fp2_ComplexPlane_multiFrekComplex.mp4'); % Save video to Desktop
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = 25; % Set the frame rate
open(v);

% Initialize time point and plot
currentDataIndex = 1; % Start at the first data point
z_offset = 0; % Offset for time fade in z-axis

% Setup time window loop
for frameNum = 1:frameCount
    clf; % Clear figure

    % Create subplot for T6
    subplot(1, 2, 1);
    hold on;

    % Plot black reference line at real=0, imag=0 for T6
    plot3([0 0], [0 0], [0 fadeLength], 'k-', 'LineWidth', 2);

    % Plot multiFrekComplex for T6
    if currentDataIndex >= fadeLength
        currentWindowT6 = T6_data(currentDataIndex-fadeLength+1:currentDataIndex);
    else
        currentWindowT6 = [zeros(fadeLength-currentDataIndex, 1); T6_data(1:currentDataIndex)];
    end
    realPartsT6 = real(currentWindowT6);
    imagPartsT6 = imag(currentWindowT6);
    z_valsT6 = z_offset + (1:length(currentWindowT6)); % Reverse the time axis
    for i = 1:length(currentWindowT6)
        opacity = (i / fadeLength);
        fadeColor = [0, 0, 1] * opacity + [1, 1, 1] * (1 - opacity); % Blue fade
        plot3([0 realPartsT6(i)], [0 imagPartsT6(i)], [z_valsT6(i) z_valsT6(i)], 'Color', fadeColor, 'LineWidth', 1);
        if i > 1
            plot3([realPartsT6(i-1) realPartsT6(i)], [imagPartsT6(i-1) imagPartsT6(i)], [z_valsT6(i-1) z_valsT6(i)], 'Color', fadeColor, 'LineWidth', 1);
        end
        plot3(realPartsT6(i), imagPartsT6(i), z_valsT6(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
    end

    xlabel('Imaginary');
    ylabel('Real');
    zlabel('Time fade');
    title('T6 Complex Plane for multiFrekComplex with Fading Vectors');
    grid on;
    view([0.1 0.1 0.2]); % Adjust this to your desired viewing angle
    axis([-20 20 -20 20 0 fadeLength]); % Reverse the time axis (time increasing upwards)

    % Create subplot for Fp2
    subplot(1, 2, 2);
    hold on;

    % Plot black reference line at real=0, imag=0 for Fp2
    plot3([0 0], [0 0], [0 fadeLength], 'k-', 'LineWidth', 2);

    % Plot multiFrekComplex for Fp2
    if currentDataIndex >= fadeLength
        currentWindowFp2 = Fp2_data(currentDataIndex-fadeLength+1:currentDataIndex);
    else
        currentWindowFp2 = [zeros(fadeLength-currentDataIndex, 1); Fp2_data(1:currentDataIndex)];
    end
    realPartsFp2 = real(currentWindowFp2);
    imagPartsFp2 = imag(currentWindowFp2);
    z_valsFp2 = z_offset + (1:length(currentWindowFp2)); % Reverse the time axis
    for i = 1:length(currentWindowFp2)
        opacity = (i / fadeLength);
        fadeColor = [0, 0, 1] * opacity + [1, 1, 1] * (1 - opacity); % Blue fade
        plot3([0 realPartsFp2(i)], [0 imagPartsFp2(i)], [z_valsFp2(i) z_valsFp2(i)], 'Color', fadeColor, 'LineWidth', 1);
        if i > 1
            plot3([realPartsFp2(i-1) realPartsFp2(i)], [imagPartsFp2(i-1) imagPartsFp2(i)], [z_valsFp2(i-1) z_valsFp2(i)], 'Color', fadeColor, 'LineWidth', 1);
        end
        plot3(realPartsFp2(i), imagPartsFp2(i), z_valsFp2(i), 'o', 'MarkerSize', 5, 'MarkerFaceColor', fadeColor, 'MarkerEdgeColor', 'none');
    end

    xlabel('Imaginary');
    ylabel('Real');
    zlabel('Time fade');
    title('Fp2 Complex Plane for multiFrekComplex with Fading Vectors');
    grid on;
    view([0.1 0.1 0.2]); % Adjust this to your desired viewing angle
    axis([-20 20 -20 20 0 fadeLength]); % Reverse the time axis (time increasing upwards)

    % Update the current data index
    currentDataIndex = currentDataIndex + stepSize;

    % Capture the frame
    frame = getframe(figh);
    writeVideo(v, frame);

    % Check if we reach the end of the data
    if currentDataIndex > size(T6_data, 1)
        break;
    end
end

% Close video writer
close(v);
close(figh); % Close the figure




%Spectral CMW buckets for each electrode table!!
% Define the electrode labels
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4'};

% List of complex number datasets for different frequencies
complexNumberSet = {complexNumbers2Hz, complexNumbers5Hz, complexNumbers7Hz, ...
                    complexNumbers10Hz, complexNumbers15Hz, complexNumbers20Hz, ...
                    complexNumbers30Hz, complexNumbers40Hz};

% Frequency labels
freqLabels = {'2 Hz', '5 Hz', '7 Hz', '10 Hz', '15 Hz', '20 Hz', '30 Hz', '40 Hz'};

% Initialize amplitude matrix
amplitudes = zeros(length(electrodeLabels), length(freqLabels));

% Calculate the amplitude for each frequency and electrode
for freqIdx = 1:length(complexNumberSet)
    currentFreqData = complexNumberSet{freqIdx};
    for ch = 1:length(electrodeLabels)
        % Calculate the amplitude (sqrt(real^2 + imag^2)) for each time point and sum
        amplitudeSum = sum(sqrt(real(currentFreqData(:, ch)).^2 + imag(currentFreqData(:, ch)).^2));
        % Store the average amplitude for the electrode and frequency
        amplitudes(ch, freqIdx) = amplitudeSum / size(currentFreqData, 1);
    end
end

% Scale the amplitude values between 0 and 1 using the maximum value across all electrodes and frequencies
maxAmplitude = max(amplitudes(:));
amplitudesScaled = amplitudes / maxAmplitude;

% Create a table with electrode labels as row names and frequency labels as column names
amplitudeTable = array2table(amplitudesScaled, 'VariableNames', freqLabels, 'RowNames', electrodeLabels);

% Display the table
disp(amplitudeTable);


% Define the flipped new electrode order (Fp1 and Fp2 at the top)
reorderedElectrodeLabelsFlipped = {'Fp1', 'Fp2', 'F7', 'F3', 'F4', 'F8', 'T3', 'C3', 'C4', 'T4', ...
                                   'T5', 'P3', 'P4', 'T6', 'O1', 'O2'};

% Find the indices of the flipped reordered labels in the original electrode labels
[~, newOrderFlipped] = ismember(reorderedElectrodeLabelsFlipped, electrodeLabels);

% Rearrange the amplitude matrix according to the flipped electrode order
amplitudesReorderedFlipped = amplitudesScaled(newOrderFlipped, :);

% Create a new table with the flipped reordered electrodes
amplitudeTableReorderedFlipped = array2table(amplitudesReorderedFlipped, 'VariableNames', freqLabels, 'RowNames', reorderedElectrodeLabelsFlipped);

% Display the flipped reordered table
disp(amplitudeTableReorderedFlipped);




% Electrode positions based on the 10-20 system (in the same order as electrodeLabels)
electrodePos = [
    -0.125, 0.5;   % Fp1
    0.125, 0.5;    % Fp2
    -0.2, 0;       % C3
    0.2, 0;        % C4
    -0.3425, -0.27;% T5 (aka P7)
    0.3425, -0.27; % T6 (aka P8)
    -0.125, -0.5;  % O1
    0.125, -0.5;   % O2
    -0.3425, 0.27; % F7
    0.3425, 0.27;  % F8
    -0.18, 0.15;   % F3
    0.18, 0.15;    % F4
    -0.5, 0;       % T3
    0.5, 0;        % T4
    -0.18, -0.15;  % P3
    0.18, -0.15    % P4
];

% Electrode labels in original order
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', 'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4'};

% Define colors for the frequencies (from purple to pink)
colors = [148 0 211; ... % Purple (2 Hz)
          75 0 130; ... % Indigo (5 Hz)
          0 0 255; ... % Blue (7 Hz)
          0 255 0; ... % Green (10 Hz)
          255 255 0; ... % Yellow (15 Hz)
          255 165 0; ... % Orange (20 Hz)
          255 0 0; ... % Red (30 Hz)
          255 105 180] / 255; % Pink (40 Hz)

% Plot the head diagram
figure;
hold on;

% Draw a rough head outline (circle for scalp)
viscircles([0 0], 1, 'LineStyle', '--', 'LineWidth', 1);

% Draw a nose representation
plot([0 0.1 -0.1 0], [1 1.1 1.1 1], 'k', 'LineWidth', 1);

% Plot the electrodes as labels on the scalp projection
for i = 1:length(electrodeLabels)
    text(electrodePos(i,1), electrodePos(i,2), electrodeLabels{i}, 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Scaling amplitude for each frequency
maxAmplitude = max(amplitudesScaled(:)); % Find the maximum amplitude for scaling

% Plot vectors for each frequency (from 2 Hz to 40 Hz)
for i = 1:length(electrodeLabels)
    for f = 1:8 % 8 frequency bands
        % Get the amplitude for the current electrode and frequency
        amplitude = amplitudesScaled(i, f) / maxAmplitude; % Normalize from 0 to 1

        % Calculate the vector direction and magnitude
        vecX = amplitude * cos(pi/4); % Direction in X (adjust if needed)
        vecY = amplitude * sin(pi/4); % Direction in Y (adjust if needed)

        % Plot the vector with the corresponding color
        quiver(electrodePos(i,1), electrodePos(i,2), vecX, vecY, 'Color', colors(f, :), 'MaxHeadSize', 0.5, 'LineWidth', 2, 'AutoScale', 'off');
    end
end

% Adjust axis for scalp projection
axis equal;
axis([-1.2 1.2 -1.2 1.3]);
set(gca, 'XTick', [], 'YTick', []);

% Add a legend
legendLabels = {'2 Hz', '5 Hz', '7 Hz', '10 Hz', '15 Hz', '20 Hz', '30 Hz', '40 Hz'};
for f = 1:8
    plot(nan, nan, 'Color', colors(f, :), 'LineWidth', 2); % Add invisible plot for legend
end
legend(legendLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');

title('Electrode Amplitude Vectors (0 to 1 magnitude)');
hold off;



%or, since not have image processing for viscircles

% Electrode positions based on the 10-20 system (in the same order as electrodeLabels)
electrodePos = [
    -0.125, 0.5;   % Fp1
    0.125, 0.5;    % Fp2
    -0.2, 0;       % C3
    0.2, 0;        % C4
    -0.3425, -0.27;% T5 (aka P7)
    0.3425, -0.27; % T6 (aka P8)
    -0.125, -0.5;  % O1
    0.125, -0.5;   % O2
    -0.3425, 0.27; % F7
    0.3425, 0.27;  % F8
    -0.18, 0.15;   % F3
    0.18, 0.15;    % F4
    -0.5, 0;       % T3
    0.5, 0;        % T4
    -0.18, -0.15;  % P3
    0.18, -0.15    % P4
];

% Electrode labels in original order
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'T5', 'T6', 'O1', 'O2', 'F7', 'F8', 'F3', 'F4', 'T3', 'T4', 'P3', 'P4'};

% Define colors for the frequencies (from purple to pink)
colors = [148 0 211; ... % Purple (2 Hz)
          75 0 130; ... % Indigo (5 Hz)
          0 0 255; ... % Blue (7 Hz)
          0 255 0; ... % Green (10 Hz)
          255 255 0; ... % Yellow (15 Hz)
          255 165 0; ... % Orange (20 Hz)
          255 0 0; ... % Red (30 Hz)
          255 105 180] / 255; % Pink (40 Hz)

% Define unique angles for each frequency (in degrees)
angles = [0, 45, 90, 135, 180, 225, 270, 315];

% Plot the head diagram
figure;
hold on;

% Draw a circle representing the head (2x bigger)
theta = linspace(0, 2*pi, 100); % Angle for circle
x_circle = cos(theta); % Scaled up to 2x size
y_circle = sin(theta);
plot(x_circle, y_circle, 'k--', 'LineWidth', 1); % Circle outline

% Draw a nose representation
plot([0 0.1 -0.1 0]*2, [1 1.1 1.1 1]*2, 'k', 'LineWidth', 1); % Scaled up nose

% Adjusted positions for Fp1 and Fp2 labels (shifted down 0.1 units)
for i = 1:length(electrodeLabels)
    if strcmp(electrodeLabels{i}, 'Fp1') || strcmp(electrodeLabels{i}, 'Fp2')
        text(electrodePos(i,1)*2, (electrodePos(i,2)*2 - 0.12), electrodeLabels{i}, 'HorizontalAlignment', 'center', 'FontSize', 12); % Shifted down
    else
        text(electrodePos(i,1)*2, (electrodePos(i,2)*2 - 0.12), electrodeLabels{i}, 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
end

% Scaling amplitude for each frequency
maxAmplitude = max(amplitudesScaled(:)); % Find the maximum amplitude for scaling

% Initialize plot handles for legend
quiverHandles = gobjects(1, 8); % Initialize handles for each frequency

% Plot vectors for each frequency (from 2 Hz to 40 Hz) with unique angles
for f = 1:8 % 8 frequency bands
    for i = 1:length(electrodeLabels)
        % Get the amplitude for the current electrode and frequency
        amplitude = amplitudesScaled(i, f) / maxAmplitude; % Normalize from 0 to 1

        % Calculate the vector direction based on the unique angle
        angleRad = deg2rad(angles(f)); % Convert to radians
        vecX = (amplitude * cos(angleRad)) / 4; % Direction in X (4x shorter)
        vecY = (amplitude * sin(angleRad)) / 4; % Direction in Y (4x shorter)

        % Plot the vector with the corresponding color
        quiverHandles(f) = quiver(electrodePos(i,1)*2, electrodePos(i,2)*2, vecX, vecY, 'Color', colors(f, :), 'MaxHeadSize', 0.5, 'LineWidth', 2, 'AutoScale', 'off');
    end
end

% Adjust axis for scalp projection
axis equal;
axis([-2.4 2.4 -2.4 2.4]); % Scaled for 2x headplot size
set(gca, 'XTick', [], 'YTick', []);

% Add a legend using the plot handles
legend(quiverHandles, {'2 Hz', '5 Hz', '7 Hz', '10 Hz', '15 Hz', '20 Hz', '30 Hz', '40 Hz'}, 'Location', 'southoutside', 'Orientation', 'horizontal');

title('Electrode Amplitude Vectors (0 to 1 magnitude)');
hold off;
