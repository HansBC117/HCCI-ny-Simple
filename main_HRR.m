clc; clear;

% === Configuration ===
filename = 'Group1Lambda6.xls';  % <-- change to other groups later
[~, sheetNames] = xlsfinfo(filename);

encoder_resolution = 0.1;  % CAD resolution
TDC_shift = 1;             % degrees
atm_pressure = 1.013;      % bar
gain = 20;                 % bar/V
gamma = 1.35;

% === Preallocate ===
HRR_matrix = zeros(length(sheetNames), 7200);  % 7200 points (0.1 CAD over 720°)
injection_timings = zeros(length(sheetNames), 1);

% === Loop over all injection sheets ===
for i = 1:length(sheetNames)
    sheet = sheetNames{i};
    data = readtable(filename, 'Sheet', sheet);

    % Extract time, fired, and motored signals
    time = data{:,1};
    fired_signal = data{:,2};
    motored_signal = data{:,3};

    % Store injection timing from sheet name (e.g., '-40', '-20')
    injection_timings(i) = str2double(sheet);

    % === Ensure signal has even number of elements for filtering ===
    if mod(length(fired_signal), 2) ~= 0
        fired_signal = fired_signal(1:end-1);
        motored_signal = motored_signal(1:end-1);
        time = time(1:end-1);
    end

    % === Filter Signals ===
    filtered_fired   = fftBPfilter(time, fired_signal,   [0 4000 1 0], 2000, 'plot_off');
    filtered_motored = fftBPfilter(time, motored_signal, [0 4000 1 0], 2000, 'plot_off');


    % === Trim signals to an integer number of engine cycles ===
    points_per_cycle = 720 / encoder_resolution;
    min_len = min(length(filtered_fired), length(filtered_motored));
    usable_len = floor(min_len / points_per_cycle) * points_per_cycle;

    filtered_fired = filtered_fired(1:usable_len);
    filtered_motored = filtered_motored(1:usable_len);
    time = time(1:usable_len);

    % === Process Pressure ===
    [~, ~, P_rep_fired, Theta, V_theta, P_fired] = ...
        process_cylinder_pressure_B(filtered_motored, filtered_fired, encoder_resolution, TDC_shift, atm_pressure, gain);
    
    [~, ~, ~, ~, ~, P_motored] = ...
        process_cylinder_pressure_B(filtered_motored, filtered_motored, encoder_resolution, TDC_shift, atm_pressure, gain);

    % === Compute HRR ===
    dTheta = mean(diff(Theta));
    dP_fired = gradient(P_fired, dTheta);
    dP_motored = gradient(P_motored, dTheta);
    dV = gradient(V_theta, dTheta);

    HRR_fired = (gamma/(gamma-1)) * P_fired .* dV + (1/(gamma-1)) * V_theta .* dP_fired;
    HRR_motored = (gamma/(gamma-1)) * P_motored .* dV + (1/(gamma-1)) * V_theta .* dP_motored;

    % === Relative HRR ===
    relative_HRR = HRR_fired - HRR_motored;

    % === Store in HRR Matrix ===
    HRR_matrix(i, :) = relative_HRR;
end

% === Plot Surface ===
CAD_vector = Theta;
[X, Y] = meshgrid(CAD_vector, injection_timings);
[Xi, Yi] = meshgrid(linspace(min(CAD_vector), max(CAD_vector), 100), ...
                    linspace(min(injection_timings), max(injection_timings), 100));
Zi = interp2(X, Y, HRR_matrix, Xi, Yi, 'spline');

figure
surf(Xi, Yi, Zi)
view([0 0 1])
shading flat
colorbar
xlabel('CAD [°]')
ylabel('Injection timing [CAD BTDC]')
zlabel('Relative HRR')
title(['Relative HRR Map - ', filename], 'Interpreter', 'none');



