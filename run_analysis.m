% Load pressure data
data = load('Exercise_A_raw_pressures.mat');

% Define input parameters for signal processing
encoder_resolution = 0.5; % deg
TDC_shift = 0.5;          % deg
atm_pressure = 1.013;     % bar
gain = 10;                % bar/V

% Call Step A processing function
[~, avg_pressure, ~, Theta, ~, ~] = process_cylinder_pressure_B( ...
    data.motored_pressure_signal, ...
    data.fired_pressure_signal, ...
    encoder_resolution, ...
    TDC_shift, ...
    atm_pressure, ...
    gain);

% Define engine and gas properties for Exercise B
S = 114.3e-3;    % stroke [m]
B = 82.6e-3;     % bore [m]
L = 254e-3;      % connecting rod [m]
CR = 7;          % compression ratio
Tinlet = 300;    % inlet temperature [K]
gamma = 1.3;     % heat capacity ratio

% Call analysis function for Programming Exercise B
[IMEP, PMEP, V_theta, AHRR, AHR, Tgas_mean] = ...
    analyze_cylinder_pressure(Theta, avg_pressure, S, B, L, CR, Tinlet, gamma);

% ---- PLOTS ----

% Plot AHRR
figure;
plot(Theta, AHRR);
xlabel('Crank Angle Degree (CAD)');
ylabel('AHRR [J/deg]');
title('Apparent Heat Release Rate');
grid on;

% Plot Cumulative Heat Release
figure;
plot(Theta, AHR);
xlabel('Crank Angle Degree (CAD)');
ylabel('Cumulative Heat Release [J]');
title('Apparent Heat Release');
grid on;

% Plot Estimated In-Cylinder Gas Temperature
figure;
plot(Theta, Tgas_mean);
xlabel('Crank Angle Degree (CAD)');
ylabel('Gas Temperature [K]');
title('Estimated Mean In-Cylinder Gas Temperature');
grid on;

% Display IMEP and PMEP
fprintf('IMEP: %.2f bar\n', IMEP / 1e5);
fprintf('PMEP: %.2f bar\n', PMEP / 1e5);