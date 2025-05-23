clc; clear;

% --- Configuration ---
% Each row: {filename, lambda value, engine speed [rpm]}
files = {
    'Group1Lambda6.xls', 6, 1200;
    'Group2Lambda6.xls', 6, 900;
    'Group3Lambda6.xlsx', 6, 1700;
    'Group4Lambda4.xls', 4, 1200;
    'Group5Lambda4.xls', 4, 900;
    'Group6Lambda4.xls', 4, 1700;
};

encoder_resolution = 0.1;  % CAD resolution
TDC_shift = 1;             % degrees
atm_pressure = 1.013;      % bar
gain = 20;                 % bar/V
gamma = 1.35;

for fIdx = 1:size(files,1)
    filename = files{fIdx,1};
    lambda   = files{fIdx,2};
    speed    = files{fIdx,3};
    [~, sheetNames] = xlsfinfo(filename);

    HRR_matrix = zeros(length(sheetNames), 7200);  % 720 deg at 0.1 CAD
    injection_timings = zeros(length(sheetNames),1);

    for i = 1:length(sheetNames)
        sheet = sheetNames{i};
        data = readtable(filename, 'Sheet', sheet);
        time = data{:,1};
        fired_signal = data{:,2};
        motored_signal = data{:,3};

        injection_timings(i) = str2double(sheet);

        % Ensure even number of samples for FFT based filtering
        if mod(length(fired_signal),2) ~= 0
            fired_signal = fired_signal(1:end-1);
            motored_signal = motored_signal(1:end-1);
            time = time(1:end-1);
        end

        % Filter signals
        filtered_fired   = fftBPfilter(time, fired_signal,   [0 4000 1 0], 2000, 'plot_off');
        filtered_motored = fftBPfilter(time, motored_signal, [0 4000 1 0], 2000, 'plot_off');


        % Trim signals to an integer number of cycles based on available data
        points_per_cycle = 720 / encoder_resolution;
        min_len = min(length(filtered_fired), length(filtered_motored));
        usable_len = floor(min_len / points_per_cycle) * points_per_cycle;

        filtered_fired = filtered_fired(1:usable_len);
        filtered_motored = filtered_motored(1:usable_len);

        % Process pressure
        [~, ~, ~, Theta, V_theta, P_fired] = ...
            process_cylinder_pressure_B(filtered_motored, filtered_fired, ...
                encoder_resolution, TDC_shift, atm_pressure, gain);
        [~, ~, ~, ~, ~, P_motored] = ...
            process_cylinder_pressure_B(filtered_motored, filtered_motored, ...
                encoder_resolution, TDC_shift, atm_pressure, gain);

        % Calculate HRR
        dTheta = mean(diff(Theta));
        dP_fired = gradient(P_fired, dTheta);
        dP_motored = gradient(P_motored, dTheta);
        dV = gradient(V_theta, dTheta);
        HRR_fired = (gamma/(gamma-1))*P_fired.*dV + ...
                    (1/(gamma-1))*V_theta.*dP_fired;
        HRR_motored = (gamma/(gamma-1))*P_motored.*dV + ...
                      (1/(gamma-1))*V_theta.*dP_motored;
        relative_HRR = HRR_fired - HRR_motored;

        HRR_matrix(i,:) = relative_HRR;
    end

    % Surface plot
    CAD_vector = Theta;
    [X,Y] = meshgrid(CAD_vector, injection_timings);
    [Xi,Yi] = meshgrid(linspace(min(CAD_vector), max(CAD_vector), 100), ...
                       linspace(min(injection_timings), max(injection_timings), 100));
    Zi = interp2(X, Y, HRR_matrix, Xi, Yi, 'spline');

    figure;
    surf(Xi, Yi, Zi);
    view([0 0 1]);
    shading flat;
    colorbar;
    xlabel('CAD [\circ]');
    ylabel('Injection timing [CAD BTDC]');
    zlabel('Relative HRR');
    title(sprintf('Relative HRR Map - \\lambda=%d, %d rpm', lambda, speed));
end