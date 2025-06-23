%% Clean and upsample OpenBCI's crappy 5-11 ms spread b/w datapoints using their at-least-accurate-timestamps and interpolation of a curve.
% Basic workflow to take voltage EEG and obtain an array
% of complex morlet wavelet coefficients for freqs 2/7/10/15/20/30/40hz
% With this array[freqs] of complex numbers, with just any
% timepoint, the complex number, you can obtain
% REAL real(eegConvress(freq,:,channel)
% IMAG imag(eegConvress(freq,:,channel)
% AMP abs(eegConvress(freq,:,channel)
% POWER abs(eegConvress(freq,:,channel).^2
% PHASE angle(eegConvress(freq,:,channel)
tableData = OpenBCIRAW20250601150442;
tableData = OpenBCIRAW20240509134543plusVapeNow;
tableData = tableData(6:end, :); %Removes the column header crap
dataArray = table2array(tableData);
interpolatedVoltages = csvread('interpolatedVoltages.csv');
electrodeLabels = {'Fp1', 'Fp2', 'C3', 'C4', 'Pz', 'Fz', 'O1', 'O2', ...
                   'F7', 'F8', 'F3', 'F4', 'T3', 'Cz', 'P3', 'P4'}; % Probably did not have Pz, Fz, or Cz on crispie test


% Given dataArray (448209 x 30 double array for large EEG data)
% Columns 1-16: Voltages, Column 30: Timepoints (unix) (Columns 17-29: crap, but had to import to get to timepoints column)

% Extract voltages (448209 x 16) and timepoints (448209 x 1)
voltages = dataArray(:, 1:16);
timepoints = dataArray(:, 30);

% Convert timepoints to milliseconds
timepoints_ms = (timepoints - min(timepoints)) * 1000;

% Define the desired sampling rate for extrapolated values
desired_interval_ms = 1; % 1 ms interval

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
        t_interpolated = unique_timepoints(1):desired_interval_ms:unique_timepoints(end);
        
        % Query the smooth curve for the extrapolated values
        voltage_interpolated = interp1(t_smooth, voltage_smooth, t_interpolated, 'spline');
        
        % Append the extrapolated values to the result for this channel
        interpolatedChannel = [interpolatedChannel; voltage_interpolated(:)];
    end
    
    % Store the result for this channel
    interpolatedVoltages(1:length(interpolatedChannel), channel) = interpolatedChannel;
end

% Ensure that interpolatedVoltages is trimmed to the correct number of points
interpolatedVoltages = interpolatedVoltages(1:num_points, :);

% Now interpolatedVoltages contains the extrapolated values for all channels


%% Sanity PLOT: check a few seconds of data to see exactly what the interpolation ended up looking like against the original
global_t = (0:(num_points - 1))' + min(timepoints_ms);
% Plot 10 seconds of data
t0 = 10 * 1000;  % ms

% Find closest indices in each dataset
[~, idx_orig] = min(abs(timepoints_ms - t0));
[~, idx_interp] = min(abs(global_t - t0));

% Extract some reasonable window around that point
window_orig = (idx_orig-50):(idx_orig+50);
window_interp = (idx_interp-400):(idx_interp+400);  % 1000 Hz vs 125 Hz

% Clean boundaries (in case windows go out of range)
window_orig = window_orig(window_orig > 0 & window_orig <= length(timepoints_ms));
window_interp = window_interp(window_interp > 0 & window_interp <= length(global_t));

% Plot
figure; hold on;

% Plot original raw data as red dots
plot(timepoints_ms(window_orig), voltages(window_orig, 1), 'ro', 'DisplayName','Original Data');

% Plot interpolated smooth curve as blue line
plot(global_t(window_interp), interpolatedVoltages(window_interp, 1), 'b-', 'LineWidth', 1.5, 'DisplayName','Interpolated Curve');

% Plot interpolated individual dots as black triangles
plot(global_t(window_interp), interpolatedVoltages(window_interp, 1), 'k^', 'MarkerSize', 4, 'DisplayName','Interpolated Dots');

xlabel('Time (ms)');
ylabel('Voltage');
title('Original vs Interpolated EEG');
legend;
grid on;
hold off;
