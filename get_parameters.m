function p = get_parameters()
% GET_PARAMETERS  Physical parameters for 3D composite box beam optimization
%
% Spar modelled as a thin-walled rectangular box beam (Euler-Bernoulli
% bending + Bredt-Batho torsion) with composite laminate walls.
%
% Design variables: x = [h, b, t, theta]
%   h     = spar box height      (mm)   [25, 70]
%   b     = spar box width       (mm)   [15, 60]
%   t     = wall thickness       (mm)   [0.5, 2.5]
%   theta = off-axis ply angle   (deg)  [0, 90]

    % --- Geometry and loading ---
    p.L        = 900;          % span length (mm)
    p.P        = 1011;         % tip vertical load (N)
    p.T_torque = 63e3;         % tip torsional moment (N*mm)  [80 Nm]

    % --- Laminate definition: [0/theta/-theta/0]s ---
    p.n_plies = 8;             % total number of plies
    p.t_ply   = 0.125;         % single ply thickness (mm)

    % --- CFRP material (Toray T300/5208) ---
    p.rho  = 1600e-9;          % density (kg/mm^3) — 1600 kg/m^3
    p.E1   = 181e3;            % longitudinal modulus     (MPa)
    p.E2   = 10.3e3;           % transverse modulus       (MPa)
    p.G12  = 7.17e3;           % in-plane shear modulus   (MPa)
    p.nu12 = 0.28;             % major Poisson's ratio

    % --- Structural allowables ---
    p.sigma_allow = 331;       % max bending stress   (MPa)
    p.tau_allow   = 70;        % max shear stress     (MPa)
    p.delta_max   = 15;        % max tip deflection   (mm)
    p.phi_max     = 2*pi/180;  % max twist angle      (rad) — 2 degrees

    % --- Design variable bounds ---
    p.lb = [25,  15,  0.5,   0];   % [h_min, b_min, t_min, theta_min]
    p.ub = [70,  60,  2.5,  90];   % [h_max, b_max, t_max, theta_max]

    % --- Variable names (for plotting) ---
    p.var_names  = {'h (mm)', 'b (mm)', 't (mm)', '\theta (deg)'};
    p.var_labels = {'h', 'b', 't', 'theta'};

end
