function smoothed_emg = smoothEMG(emg_data, window_size)
    % Check if the input is a valid vector
    if ~isvector(emg_data)
        error('Input must be a vector.');
    end
    
    % Check if the window size is valid
    if window_size < 1 || window_size > length(emg_data)
        error('Invalid window size.');
    end
    
    % Calculate the RMS over a moving window
    half_window = floor(window_size / 2);
    smoothed_emg = zeros(size(emg_data));
    
    for i = 1:length(emg_data)
        start_idx = max(1, i - half_window);
        end_idx = min(length(emg_data), i + half_window);
        window_data = emg_data(start_idx:end_idx);
        smoothed_emg(i) = sqrt(mean(window_data .^ 2));
    end
end
