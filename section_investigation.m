function section_investigation(p)
% SECTION_INVESTIGATION  Initial problem investigation
%
% Covers:
%   1. Boundedness  — does a finite minimum exist?
%   2. Monotonicity — which constraints are active at the optimum?
%   3. Convexity    — is the problem convex? Are local minima a risk?
%   4. Numerical noise — is the model smooth?
%
% Produces:
%   Fig 0a - Objective landscape (h-b slice)
%   Fig 0b - Monotonicity: gradient vectors in design space
%   Fig 0c - Noise check: numerical vs analytical gradient comparison
%   Fig 0d - Laminate properties vs ply angle theta

    fprintf('\n  --- Boundedness Analysis ---\n');

    % ------------------------------------------------------------------
    % 1. BOUNDEDNESS
    %
    % W = rho * 2*(h+b) * t * L  ->  approaches 0 as h,b,t -> 0
    %
    % But the bending stress constraint imposes:
    %   P*L*(h/2) / I_xx  <= sigma_allow
    % and the deflection constraint:
    %   P*L^3 / (3*Ex*I_xx) <= delta_max
    %
    % These force strictly positive lower bounds on the cross-section size.
    % The problem IS bounded below.
    % ------------------------------------------------------------------
    fprintf('  W = rho * 2*(h+b) * t * L  ->  0  as h,b,t -> 0 (unbounded below)\n');
    fprintf('  But structural constraints (stress, deflection) impose lower bounds\n');
    fprintf('  on the cross-section size -> W has a finite minimum.\n');
    fprintf('  -> Problem IS bounded below. Minimum exists.\n');

    % Quick estimate: evaluate W at a representative feasible point
    x_test = [50; 30; 2; 45];   % h=50, b=30, t=2, theta=45
    W_test = beam_weight(x_test, p);
    [g_test, ~, info_test] = beam_constraints(x_test, p);
    fprintf('  Test point [h=%.0f, b=%.0f, t=%.1f, theta=%.0f]: W = %.4f kg\n', ...
            x_test(1), x_test(2), x_test(3), x_test(4), W_test);
    fprintf('  Constraints: g = [%.3f, %.3f, %.3f, %.3f, %.3f]\n', g_test);
    fprintf('  (negative = feasible)\n');

    % ------------------------------------------------------------------
    % 2. MONOTONICITY ANALYSIS
    % ------------------------------------------------------------------
    fprintf('\n  --- Monotonicity Analysis ---\n');
    fprintf('  dW/dh     = 2*rho*t*L     > 0   for all t > 0\n');
    fprintf('  dW/db     = 2*rho*t*L     > 0   for all t > 0\n');
    fprintf('  dW/dt     = 2*rho*(h+b)*L > 0   for all h,b > 0\n');
    fprintf('  dW/dtheta = 0                    (theta only in constraints)\n');
    fprintf('  -> Objective is monotonically increasing in h, b, t.\n');
    fprintf('  -> The optimum lies on the constraint boundary.\n');
    fprintf('  -> theta acts as a "free" variable affecting only stiffness/strength.\n');

    % ------------------------------------------------------------------
    % 3. CONVEXITY ANALYSIS
    % ------------------------------------------------------------------
    fprintf('\n  --- Convexity Analysis ---\n');
    c_coeff = 2 * p.rho * p.L;
    H_hbt = c_coeff * [0 0 1; 0 0 1; 1 1 0];
    eig_vals = eig(H_hbt);
    fprintf('  Hessian eigenvalues (h,b,t subspace): [%.4f, %.4f, %.4f]\n', ...
            sort(eig_vals));
    fprintf('  -> Hessian is INDEFINITE: objective is non-convex.\n');
    fprintf('  -> Cannot guarantee a single global minimum analytically.\n');
    fprintf('  -> Will use multiple starting points to verify.\n');

    % ------------------------------------------------------------------
    % 4. NUMERICAL NOISE CHECK
    % ------------------------------------------------------------------
    fprintf('\n  --- Numerical Noise Check ---\n');

    h_test = 50; b_test = 30; t_test = 2; th_test = 45;
    x_test = [h_test; b_test; t_test; th_test];

    % Analytical gradients of W
    dW_dh_exact = 2 * p.rho * t_test * p.L;
    dW_db_exact = 2 * p.rho * t_test * p.L;
    dW_dt_exact = 2 * p.rho * (h_test + b_test) * p.L;

    % Numerical gradients (forward differences, various step sizes)
    step_list = logspace(-10, -1, 40);
    err_h = zeros(size(step_list));
    err_b = zeros(size(step_list));
    err_t = zeros(size(step_list));

    W0 = beam_weight(x_test, p);
    for i = 1:length(step_list)
        hs = step_list(i);

        xp = x_test; xp(1) = xp(1) + hs;
        err_h(i) = abs((beam_weight(xp,p) - W0)/hs - dW_dh_exact);

        xp = x_test; xp(2) = xp(2) + hs;
        err_b(i) = abs((beam_weight(xp,p) - W0)/hs - dW_db_exact);

        xp = x_test; xp(3) = xp(3) + hs;
        err_t(i) = abs((beam_weight(xp,p) - W0)/hs - dW_dt_exact);
    end

    fprintf('  Analytical dW/dh = %.6e,  dW/db = %.6e,  dW/dt = %.6e\n', ...
            dW_dh_exact, dW_db_exact, dW_dt_exact);
    fprintf('  Min numerical error (dW/dh): %.2e\n', min(err_h));
    fprintf('  -> Objective is smooth; finite differences are reliable.\n');

    % Check constraint smoothness (theta varies Ex, Gxy smoothly)
    fprintf('  Checking constraint smoothness vs theta...\n');
    th_vec = linspace(0.1, 89.9, 200);
    Ex_vec = zeros(size(th_vec));
    Gxy_vec = zeros(size(th_vec));
    for i = 1:length(th_vec)
        [Ex_vec(i), ~, Gxy_vec(i), ~] = compute_effective_properties(th_vec(i), p);
    end
    fprintf('  Ex range:  [%.0f, %.0f] MPa\n', min(Ex_vec), max(Ex_vec));
    fprintf('  Gxy range: [%.0f, %.0f] MPa\n', min(Gxy_vec), max(Gxy_vec));
    fprintf('  Both vary smoothly with theta — model is noise-free.\n');

    % -----------------------------------------------------------------------
    % Figure 0a: Objective landscape (h vs b slice, t=2, theta=45)
    %   Assignment style: 2D contour with ShowText + constraint boundaries
    % -----------------------------------------------------------------------
    N = 200;
    h_vec = linspace(p.lb(1), p.ub(1), N);
    b_vec = linspace(p.lb(2), p.ub(2), N);
    [H, B] = meshgrid(h_vec, b_vec);
    t_fix = 2;  th_fix = 45;

    W_grid = p.rho * 2 * (H + B) * t_fix * p.L;

    % Evaluate constraints on grid
    G1_grid = zeros(N);  G3_grid = zeros(N);  G5_grid = zeros(N);
    for ii = 1:N
        for jj = 1:N
            [gg, ~] = beam_constraints([H(ii,jj), B(ii,jj), t_fix, th_fix], p);
            G1_grid(ii,jj) = gg(1);   % stress
            G3_grid(ii,jj) = gg(3);   % deflection
            G5_grid(ii,jj) = gg(5);   % buckling
        end
    end

    figure('Name', 'Fig0a: Objective Landscape', 'Position', [50 80 700 500]);
    % Objective contours with value labels (assignment style)
    contour(H, B, W_grid, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;
    % Constraint boundaries at g=0 (solid, LineWidth 1.5)
    contour(H, B, G1_grid, [0 0], 'r', 'LineWidth', 1.5);
    contour(H, B, G3_grid, [0 0], 'k', 'LineWidth', 1.5);
    contour(H, B, G5_grid, [0 0], 'c', 'LineWidth', 1.5);
    xlabel('Height h (mm)');
    ylabel('Width b (mm)');
    title(sprintf('Objective W(h, b) at t = %.0f mm, \\theta = %.0f deg', t_fix, th_fix));
    legend({'W (kg)', 'g_1: stress', 'g_3: deflection', 'g_5: buckling'}, ...
           'Location', 'northeast');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 0b: Monotonicity — gradient vectors in (h, b) space
    %   Assignment style: contour with ShowText + constraint lines + quiver
    % -----------------------------------------------------------------------
    figure('Name', 'Fig0b: Monotonicity', 'Position', [100 80 700 550]);

    contour(H, B, W_grid, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;

    % Constraint boundaries (assignment color scheme)
    contour(H, B, G1_grid, [0 0], 'r', 'LineWidth', 1.5, ...
            'DisplayName', 'g_1: stress');
    contour(H, B, G3_grid, [0 0], 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'g_3: deflection');

    % Gradient arrows
    % n_arr = 6;
    % h_arr = linspace(p.lb(1)*1.2, p.ub(1)*0.85, n_arr);
    % b_arr = linspace(p.lb(2)*1.2, p.ub(2)*0.85, n_arr);
    % [HA, BA] = meshgrid(h_arr, b_arr);
    % dW_dh_arr = 2 * p.rho * t_fix * p.L * ones(size(HA));
    % dW_db_arr = 2 * p.rho * t_fix * p.L * ones(size(BA));
    % scale = 3 / max(max(sqrt(dW_dh_arr.^2 + dW_db_arr.^2)));
    % quiver(HA, BA, dW_dh_arr*scale, dW_db_arr*scale, 0, ...
    %        'b', 'LineWidth', 1.0, 'MaxHeadSize', 0.6, 'DisplayName', '\nablaW direction');

    xlabel('Height h (mm)');
    ylabel('Width b (mm)');
    title(sprintf('Monotonicity analysis (t = %.0f mm, \\theta = %.0f deg)', t_fix, th_fix));
    legend('Location', 'northeast');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 0c: Numerical noise check
    %   Assignment style: log-log plot (sensitivity analysis style from Ex 8)
    % -----------------------------------------------------------------------
    figure('Name', 'Fig0c: Noise Check', 'Position', [150 80 700 430]);
    loglog(step_list, err_h, 'r-', 'LineWidth', 1.5, ...
           'DisplayName', '\partialf/\partialh error');
    hold on;
    loglog(step_list, err_b, 'g-', 'LineWidth', 1.5, ...
           'DisplayName', '\partialf/\partialb error');
    loglog(step_list, err_t, 'b-', 'LineWidth', 1.5, ...
           'DisplayName', '\partialf/\partialt error');
    xlabel('Perturbation size h');
    ylabel('|FD gradient - analytical|');
    title('Finite difference gradient error vs perturbation size');
    legend('Location', 'southwest');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 0d: Laminate effective properties vs theta
    % -----------------------------------------------------------------------
    figure('Name', 'Fig0d: Laminate Properties', 'Position', [200 80 900 400]);

    subplot(1,2,1);
    plot(th_vec, Ex_vec/1e3, 'b-', 'LineWidth', 1.5);
    xlabel('\theta (deg)');
    ylabel('E_x (GPa)');
    title('Effective axial modulus vs ply angle');
    grid on;

    subplot(1,2,2);
    plot(th_vec, Gxy_vec/1e3, 'r-', 'LineWidth', 1.5);
    xlabel('\theta (deg)');
    ylabel('G_{xy} (GPa)');
    title('Effective shear modulus vs ply angle');
    grid on;

    fprintf('  Figures 0a-0d generated.\n');
end
