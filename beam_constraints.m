function [c, ceq, info] = beam_constraints(x, p)
% BEAM_CONSTRAINTS  Nonlinear inequality constraints for the 3D composite box beam
%
%   Returns c <= 0 (feasible when negative)
%
%   g(1): Bending stress       sigma_bend / sigma_allow - 1
%   g(2): Torsional shear      tau_torsion / tau_allow  - 1
%   g(3): Tip deflection       delta / delta_max        - 1
%   g(4): Twist angle          phi / phi_max            - 1
%   g(5): Web panel buckling   sigma_bend / sigma_cr    - 1
%
%   x = [h, b, t, theta]

    h     = x(1);
    b     = x(2);
    t     = x(3);
    theta = x(4);

    % --- Effective laminate properties (CLT) ---
    [Ex, ~, Gxy, nu_xy] = compute_effective_properties(theta, p);

    % --- Cross-section properties (rectangular hollow section) ---
    h_in = h - 2*t;   % inner height
    b_in = b - 2*t;   % inner width

    % Guard against invalid geometry (t too large)
    if h_in <= 0 || b_in <= 0
        c    = ones(1, 5) * 1e6;   % heavily infeasible
        ceq  = [];
        if nargout > 2, info = struct(); end
        return;
    end

    I_xx  = (b*h^3 - b_in*h_in^3) / 12;        % second moment of area (mm^4)

    % Enclosed area and perimeter along mid-line (Bredt-Batho)
    h_mid = h - t;
    b_mid = b - t;
    A_enc = h_mid * b_mid;                        % enclosed area (mm^2)
    J     = 2 * A_enc^2 * t / (h_mid + b_mid);   % torsion constant (mm^4)

    % --- Internal loads at root ---
    M = p.P * p.L;              % bending moment  (N*mm)
    T = p.T_torque;             % torque           (N*mm)

    % --- Stress responses ---
    sigma_bend  = M * (h/2) / I_xx;              % max bending stress (MPa)
    tau_torsion = T / (2 * A_enc * t);           % shear stress       (MPa)

    % --- Displacement responses ---
    delta = p.P * p.L^3 / (3 * Ex * I_xx);      % tip deflection (mm)
    phi   = T * p.L     / (Gxy * J);            % twist angle    (rad)

    % --- Local buckling of web panel (simply-supported plate) ---
    k_buck   = 4.0;   % buckling coefficient
    sigma_cr = k_buck * pi^2 * Ex * (t / h_in)^2 / (12 * (1 - nu_xy^2));

    % --- Normalised constraints: g_i <= 0 is feasible ---
    c(1) = sigma_bend  / p.sigma_allow - 1;   % bending stress
    c(2) = tau_torsion / p.tau_allow   - 1;   % torsional shear
    c(3) = delta       / p.delta_max   - 1;   % tip deflection
    c(4) = phi         / p.phi_max     - 1;   % twist angle
    c(5) = sigma_bend  / sigma_cr      - 1;   % web buckling

    ceq = [];   % no equality constraints

    % --- Detailed info for post-processing ---
    if nargout > 2
        info.sigma_bend  = sigma_bend;
        info.tau_torsion = tau_torsion;
        info.delta       = delta;
        info.phi_deg     = phi * 180/pi;
        info.sigma_cr    = sigma_cr;
        info.Ex          = Ex;
        info.Gxy         = Gxy;
        info.nu_xy       = nu_xy;
        info.I_xx        = I_xx;
        info.J           = J;
    end

end
