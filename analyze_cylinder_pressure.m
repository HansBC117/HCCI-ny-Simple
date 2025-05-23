function [IMEP, PMEP, V_theta, AHRR, AHR, T_gas] = ...
    analyze_cylinder_pressure(Theta, avg_pressure_cycle, S, B, L, CR, T_inlet, gamma)

    R = S / 2;
    V_swept = (pi * B^2 / 4) * S;
    V_c = V_swept / (CR - 1);

    theta_rad = deg2rad(Theta);
    x = R * (1 - cos(theta_rad)) + ...
        (R^2 / L) * (1 - sqrt(1 - ((L/R)^-2) * sin(theta_rad).^2));
    V_theta = V_c + (pi * B^2 / 4) .* x;

    P = avg_pressure_cycle * 1e5;  % Pa
    dV = gradient(V_theta);
    W = trapz(P .* dV);
    IMEP = W / V_swept;

    intake_mask = Theta < -180 | Theta > 180;
    PMEP = trapz(P(intake_mask) .* dV(intake_mask)) / V_swept;

    dP = gradient(P);
    dTheta = gradient(theta_rad);
    AHRR = (gamma / (gamma - 1)) * P .* (dV ./ dTheta) + ...
           (1 / (gamma - 1)) * V_theta .* (dP ./ dTheta);
    AHR = cumtrapz(theta_rad, AHRR);

    R_gas = 287;  % J/kg/K
    T_gas = (P .* V_theta) / R_gas;
end
