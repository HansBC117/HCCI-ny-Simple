function [corrected_cycles, avg_pressure_cycle, rep_pressure_cycle, Theta, V_theta, P_rep] = process_cylinder_pressure_B( ...
    motored_pressure_signal, fired_pressure_signal, encoder_resolution, TDC_shift, atm_pressure, gain)

    % Step 3: Determine how many complete cycles are present
    points_per_cycle = 720 / encoder_resolution; % 720 degrees per 4-stroke cycle
    min_len = min(length(motored_pressure_signal), length(fired_pressure_signal));
    num_cycles = floor(min_len / points_per_cycle);

    if num_cycles < 1
        error('Input signal does not contain a full engine cycle.');
    end

    usable_len = num_cycles * points_per_cycle;

    % Reshape the motored pressure signal into a matrix of individual cycles
    motored_cycles = reshape(motored_pressure_signal(1:usable_len), points_per_cycle, num_cycles);
    
    % Step 4: Compute the average motored pressure cycle
    avg_motored_pressure = mean(motored_cycles, 2);
    
    % Step 5: Find the index of the maximum pressure in the motored signal
    [~, max_idx] = max(avg_motored_pressure);
    
    % Step 6: Compute the TDC index by shifting based on the given TDC shift
    TDC_idx = max_idx + round(TDC_shift / encoder_resolution);
    
    % Step 7: Reshape the fired pressure signal into cycles
    fired_cycles = reshape(fired_pressure_signal(1:usable_len), points_per_cycle, num_cycles);
    
    % Step 7: Find the max pressure index in the fired signal
    [~, max_fired_idx] = max(mean(fired_cycles, 2));
    
    % Step 8: Determine the shift required to align high-pressure TDC at 0 degrees
    shift_amount = TDC_idx - max_fired_idx;
    fired_cycles = circshift(fired_cycles, shift_amount);
    
    % Step 9: Convert pressure signals from voltage to absolute pressure (bar)
    motored_cycles = motored_cycles * gain; 
    fired_cycles = fired_cycles * gain;
    
    % Step 10: Convert the fired pressure signal into absolute pressure
    % Assume that pressure at BDC (Bottom Dead Center) is atmospheric
    BDC_index = round((540 / encoder_resolution)); % BDC is at 540 degrees
    offset = atm_pressure - mean(fired_cycles(BDC_index, :));
    corrected_cycles = fired_cycles + offset;
    
    % Step 11: Compute the average pressure cycle
    avg_pressure_cycle = mean(corrected_cycles, 2);
    
    % Step 12: Compute the most representative cycle (smallest deviation from average)
    % Compute standard deviation for each cycle
    std_dev = std(corrected_cycles, 0, 2);

    % Find the most representative cycle (smallest deviation from mean)
    [~, min_std_idx] = min(mean(std_dev, 1)); % Fix: Use dimension 1 to get correct indexing

    % Ensure min_std_idx does not exceed the number of cycles
    min_std_idx = min(max(min_std_idx, 1), size(corrected_cycles, 2)); 

    % Select the representative cycle
        rep_pressure_cycle = corrected_cycles(:, min_std_idx);
    
    % Step 13: Generate the crank angle vector from -360 to 359.5 degrees
    Theta = linspace(-360, 359.5, points_per_cycle);

     % Step 14: Function is complete; it can now be called in a test script.

    % Step 15: Display key results for debugging
    fprintf('Number of corrected cycles: %d\n', size(corrected_cycles, 2));
    fprintf('Max pressure in average cycle: %.2f bar\n', max(avg_pressure_cycle));
    fprintf('Max pressure in representative cycle: %.2f bar\n', max(rep_pressure_cycle));
    fprintf('TDC shift applied: %.2f degrees\n', TDC_shift);
    fprintf('Crank angle range: %.2f to %.2f degrees\n', Theta(1), Theta(end));
    fprintf('Atmospheric pressure used for correction: %.3f bar\n', atm_pressure);
    
    
    % ----- Step 16: Volume & Pressure for p-V plots -----
    % Engine geometry
    B = 82.6e-3;     % Bore [m]
    S = 114.3e-3;    % Stroke [m]
    L = 254e-3;      % Connecting rod [m]
    R = S / 2;       % Crank radius
    r = 7;           % Compression ratio

    V_swept = (pi * (B^2) / 4) * S;
    V_c = V_swept / (r - 1);

    theta_rad = deg2rad(Theta);
    x_theta = R * (1 - cos(theta_rad)) + ...
              (R^2 / L) * (1 - sqrt(1 - ((L / R)^-2) * (sin(theta_rad)).^2));
    V_theta = V_c + (pi * B^2 / 4) .* x_theta;

    % Convert representative pressure cycle to Pascals
    P_rep = rep_pressure_cycle * 1e5;  % Pa
end









