function section4_simplified(p)
% SECTION4_SIMPLIFIED  Simplified 2-D and 1-D optimization
%
% Reduces the 4-variable problem by fixing: theta & h 
%
% Then optimizes over (b, t) in 2D, and over t alone in 1D.
% Uses fminsearch (Nelder-Mead, derivative-free) with exterior penalty.

    h_fixed  = 60;    % fixed height  (mm)
    th_fixed = 60;    % fixed theta   (deg)
    fprintf('  Fixed: h = %.0f mm, theta = %.0f deg\n', h_fixed, th_fixed);

    % =====================================================================
    % Part A: 2-D simplified (b, t)
    % =====================================================================
    fprintf('\n  --- 2-D Simplified: optimize (b, t) ---\n');
    
    % 2D functions
    obj2D = @(x2) beam_weight([h_fixed, x2(1), x2(2), th_fixed], p);
    con2D = @(x2) beam_constraints([h_fixed, x2(1), x2(2), th_fixed], p);
    
    lb2 = [p.lb(2); p.lb(3)];
    ub2 = [p.ub(2); p.ub(3)];
    
    % fmincon
    opts_fmincon = optimoptions('fmincon', 'Algorithm', 'sqp', ...
        'Display', 'off', 'OptimalityTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-12, ...
        'MaxIterations', 1000);
    
    n_starts = 8;
    rng(42);
    b0_list = lb2(1) + rand(n_starts,1) * (ub2(1) - lb2(1));
    t0_list = lb2(2) + rand(n_starts,1) * (ub2(2) - lb2(2));
    
    fprintf('  %-10s %-10s %-10s %-10s %-10s\n', 'b0','t0','b*','t*','W*');
    fprintf('  %s\n', repmat('-', 1, 54));
    
    best_W = Inf;  b_best = NaN;  t_best = NaN;
    
    % Storage for trajectory plotting
    conv_b    = NaN(n_starts,1);  conv_t    = NaN(n_starts,1);
    conv_feas = false(n_starts,1);
    
    for i = 1:n_starts
        x0 = [b0_list(i); t0_list(i)];
        try
            [x2_opt, W2_opt, ef] = fmincon(obj2D, x0, [],[],[],[], lb2, ub2, ...
                @(x) deal(con2D(x), []), opts_fmincon);
            [c2, ~] = con2D(x2_opt);
            is_feas = ef > 0 && all(c2 <= 1e-4);
    
            conv_b(i)    = x2_opt(1);
            conv_t(i)    = x2_opt(2);
            conv_feas(i) = is_feas;
    
            if is_feas && W2_opt < best_W
                best_W = W2_opt;
                b_best = x2_opt(1);
                t_best = x2_opt(2);
            end
            if is_feas
                fprintf('  %-10.2f %-10.2f %-10.2f %-10.2f %-10.6f\n', ...
                        x0(1), x0(2), x2_opt(1), x2_opt(2), W2_opt);
            end
        catch
        end
    end
    
    fprintf('\n  Best 2D result: b* = %.2f mm, t* = %.2f mm, W* = %.6f kg\n', ...
            b_best, t_best, best_W);
    
    % -----------------------------------------------------------------------
    % Figure 7: 2D design space (b vs t)
    % -----------------------------------------------------------------------
    N = 200;
    b_vec = linspace(p.lb(2), p.ub(2), N);
    t_vec = linspace(p.lb(3), p.ub(3), N);
    [Bg, Tg] = meshgrid(b_vec, t_vec);
    
    W2D = zeros(N);  G2D = zeros(N, N, 5);
    for ii = 1:N
        for jj = 1:N
            W2D(ii,jj) = beam_weight([h_fixed, Bg(ii,jj), Tg(ii,jj), th_fixed], p);
            [gg,~] = beam_constraints([h_fixed, Bg(ii,jj), Tg(ii,jj), th_fixed], p);
            G2D(ii,jj,:) = gg;
        end
    end
    
    figure('Name', 'Fig7: 2D Simplified', 'Position', [50 400 750 560]);
    contour(Bg, Tg, W2D, 15, 'ShowText', 'on', 'LabelSpacing', 500, 'DisplayName','W (kg)');
    hold on;
    
    % Constraint boundaries
    contour(Bg, Tg, G2D(:,:,1), [0 0], 'r',  'LineWidth', 1.5, 'DisplayName', 'g_1: Bending stress');
    contour(Bg, Tg, G2D(:,:,3), [0 0], 'k',  'LineWidth', 1.5, 'DisplayName', 'g_3: Tip deflection');
    contour(Bg, Tg, G2D(:,:,4), [0 0], 'b',  'LineWidth', 1.5, 'DisplayName', 'g_4: Twist angle');
    contour(Bg, Tg, G2D(:,:,5), [0 0], 'c',  'LineWidth', 1.5, 'DisplayName', 'g_5: Web buckling');
    
    % Dashed infeasible-side indicators
    contour(Bg, Tg, G2D(:,:,1), [0.05 0.05], '--r', 'LineWidth', 1, 'HandleVisibility', 'off');
    contour(Bg, Tg, G2D(:,:,3), [0.05 0.05], '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    contour(Bg, Tg, G2D(:,:,4), [0.05 0.05], '--b', 'LineWidth', 1, 'HandleVisibility', 'off');
    contour(Bg, Tg, G2D(:,:,5), [0.05 0.05], '--c', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Single starting point (use first start, or whichever you prefer)
    i_show = 1;
    plot(b0_list(i_show), t0_list(i_show), 'o', ...
         'Color', [0.2 0.2 0.8], 'MarkerFaceColor', [0.5 0.5 1.0], ...
         'MarkerSize', 9, 'DisplayName', 'Starting point');
    
    % Best optimum (filled red pentagram)
    if ~isnan(b_best)
        plot(b_best, t_best, 'rp', 'MarkerSize', 14, 'LineWidth', 1.5, ...
             'MarkerFaceColor', 'r', 'DisplayName', 'Optimum');
    end
    
    xlabel('Width b (mm)', 'FontWeight','bold');
    ylabel('Thickness t (mm)', 'FontWeight','bold');
    title(sprintf('2-D design space (h = %.0f mm, \\theta = %.0f deg)', h_fixed, th_fixed));
    legend('Location', 'northeast');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 7b: Local minima sensitivity — multiple starting points
    % -----------------------------------------------------------------------
    fprintf('\n  --- Local minima sensitivity (grid + random starts) ---\n');
    
    % --- Define starting points ---
    % Grid starts (4x4)
    n_grid  = 4;
    b_grid  = linspace(lb2(1)*1.05, ub2(1)*0.95, n_grid);
    t_grid  = linspace(lb2(2)*1.05, ub2(2)*0.95, n_grid);
    [Bg0, Tg0] = meshgrid(b_grid, t_grid);
    b_starts_grid = Bg0(:);
    t_starts_grid = Tg0(:);
    
    % Random starts (16 points, different seed)
    n_rand = 16;
    rng(99);
    b_starts_rand = lb2(1) + rand(n_rand,1) * (ub2(1) - lb2(1));
    t_starts_rand = lb2(2) + rand(n_rand,1) * (ub2(2) - lb2(2));
    
    % Combine
    b_starts = [b_starts_grid; b_starts_rand];
    t_starts = [t_starts_grid; t_starts_rand];
    n_total  = length(b_starts);
    is_grid  = [true(length(b_starts_grid),1); false(n_rand,1)];
    
    % --- Run fmincon from each starting point ---
    opt_b     = NaN(n_total,1);
    opt_t     = NaN(n_total,1);
    opt_W     = NaN(n_total,1);
    opt_feas  = false(n_total,1);
    
    fprintf('  Running fmincon from %d starting points...\n', n_total);
    for i = 1:n_total
        x0 = [b_starts(i); t_starts(i)];
        try
            [x_opt, W_opt, ef] = fmincon(obj2D, x0, [],[],[],[], lb2, ub2, ...
                @(x) deal(con2D(x), []), opts_fmincon);
            [c_opt, ~] = con2D(x_opt);
            is_feas = ef > 0 && all(c_opt <= 1e-4);
            opt_b(i)    = x_opt(1);
            opt_t(i)    = x_opt(2);
            opt_W(i)    = W_opt;
            opt_feas(i) = is_feas;
        catch
        end
    end
    
    % Report unique optima
    feas_idx = find(opt_feas);
    if ~isempty(feas_idx)
        opt_pts = round([opt_b(feas_idx), opt_t(feas_idx)], 2);
        unique_opts = unique(opt_pts, 'rows');
        fprintf('  Unique feasible optima found: %d\n', size(unique_opts,1));
        for k = 1:size(unique_opts,1)
            fprintf('    b* = %.2f mm,  t* = %.2f mm\n', unique_opts(k,1), unique_opts(k,2));
        end
    end
    
    % --- Plot ---
    figure('Name', 'Fig7b: Local Minima Sensitivity', 'Position', [100 350 800 580]);
    contour(Bg, Tg, W2D, 15, 'ShowText', 'on', 'LabelSpacing', 500, 'DisplayName','W (kg)');
    hold on;
    
    % Constraint boundaries (same as Fig 7)
    contour(Bg, Tg, G2D(:,:,1), [0 0], 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    contour(Bg, Tg, G2D(:,:,3), [0 0], 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    contour(Bg, Tg, G2D(:,:,4), [0 0], 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    contour(Bg, Tg, G2D(:,:,5), [0 0], 'c', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Grid starting points (blue squares)
    plot(b_starts(is_grid),  t_starts(is_grid),  's', ...
         'Color', [0.1 0.1 0.8], 'MarkerFaceColor', [0.5 0.5 1.0], ...
         'MarkerSize', 8, 'DisplayName', 'Start: grid');
    
    % Random starting points (blue circles)
    plot(b_starts(~is_grid), t_starts(~is_grid), 'o', ...
         'Color', [0.1 0.1 0.8], 'MarkerFaceColor', [0.6 0.8 1.0], ...
         'MarkerSize', 7, 'DisplayName', 'Start: random');
    
    % Feasible optima (green pentagrams)
    if any(opt_feas)
        plot(opt_b(opt_feas),  opt_t(opt_feas),  'p', ...
             'Color', [0.0 0.5 0.0], 'MarkerFaceColor', [0.2 0.8 0.2], ...
             'MarkerSize', 11, 'DisplayName', 'Optimum (feasible)');
    end
    
    % Infeasible optima (red crosses)
    if any(~opt_feas & ~isnan(opt_b))
        plot(opt_b(~opt_feas & ~isnan(opt_b)), ...
             opt_t(~opt_feas & ~isnan(opt_b)), 'rx', ...
             'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Optimum (infeasible)');
    end
    
    xlabel('Width b (mm)', 'FontWeight','bold');
    ylabel('Thickness t (mm)', 'FontWeight','bold');
    title(sprintf('Local minima sensitivity (h = %.0f mm, \\theta = %.0f deg)', ...
          h_fixed, th_fixed));
    legend('Location', 'northeast');
    grid on; hold off;

    % =====================================================================
    % Part B: 1-D simplified (t only)
    % =====================================================================
    b_fixed = b_best;   % use optimal b from 2D result
    if isnan(b_fixed), b_fixed = 30; end   % fallback

    fprintf('\n  --- 1-D Simplified: optimize t only ---\n');
    fprintf('  Fixed: h = %.0f mm, b = %.2f mm, theta = %.0f deg\n', ...
            h_fixed, b_fixed, th_fixed);

    obj1D   = @(t) beam_weight([h_fixed, b_fixed, t, th_fixed], p);
    g_all   = @(t) beam_constraints([h_fixed, b_fixed, t, th_fixed], p);

    % Exterior penalty (1D)
    mu = 1e7;
    t_range = p.ub(3) - p.lb(3);
    obj_pen = @(t) obj1D(t) + mu * sum(max(0, g_all(t)).^2) + ...
              mu * (max(0, (p.lb(3)-t)/t_range)^2 + max(0, (t-p.ub(3))/t_range)^2);

    % -----------------------------------------------------------------------
    % Figure 8: 1D landscape
    % -----------------------------------------------------------------------
    t_1d   = linspace(p.lb(3), p.ub(3), 600);
    W_1d   = arrayfun(obj1D, t_1d);
    G_1d   = zeros(length(t_1d), 5);
    for i = 1:length(t_1d)
        [gg, ~] = g_all(t_1d(i));
        G_1d(i,:) = gg;
    end
    feas_1d = all(G_1d <= 0, 2);

    figure('Name', 'Fig8: 1D Landscape', 'Position', [100 400 800 560]);

    subplot(2,1,1);
    t_f = t_1d(feas_1d);
    if ~isempty(t_f)
        % ylims = [0, max(W_1d)*1.15];
        patch([t_f, fliplr(t_f)], ...
              [zeros(size(t_f)), fliplr(W_1d(feas_1d))], ...
              [0.6 1.0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
              'DisplayName', 'Feasible region');
    end
    hold on;
    plot(t_1d, W_1d, 'b-', 'LineWidth', 1.5, 'DisplayName', 'W(t)');
    xlabel('Thickness t (mm)', 'FontWeight', 'bold');
    ylabel('Weight W (kg)', 'FontWeight', 'bold');
    title(sprintf('1-D objective (h=%.0f, b=%.1f, \\theta=%.0f)', ...
          h_fixed, b_fixed, th_fixed));
    legend('Location', 'northwest'); grid on; hold off;

    subplot(2,1,2);
    % Assignment color scheme for constraints
    colors = {'r-', 'g-', 'k-', 'b-', 'c-'};
    labels = {'g_1: stress', 'g_2: shear', 'g_3: deflection', ...
              'g_4: twist', 'g_5: buckling'};
    for j = 1:5
        plot(t_1d, G_1d(:,j), colors{j}, 'LineWidth', 1.5, 'DisplayName', labels{j});
        hold on;
    end
    yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'g = 0 boundary');
    xlabel('Thickness t (mm)', 'FontWeight','bold');
    ylabel('g(t)', 'FontWeight','bold');
    title('Constraint values (g \leq 0 is feasible)');
    legend('Location', 'northeast'); grid on; hold off;

    % -----------------------------------------------------------------------
    % Part C: Algorithm Comparison on 1D penalised problem
    % -----------------------------------------------------------------------
    fprintf('\n  --- Self-implemented Optimization (1-D) ---\n');
    
    % 1. Run Standard GSS
    % (Note: ensure golden_section_search.m is updated to return hist_gss)
    [t_gss, ~, n_gss, hist_gss] = golden_section_search(obj_pen, p.lb(3), p.ub(3), 1e-10);
    
    % 2. Run Custom Hybrid
    [t_hyb, ~, n_hyb, hist_hyb] = custom_hybrid_search(obj_pen, p.lb(3), p.ub(3), 1e-10);
    W_hyb = obj1D(t_hyb);
    [gg_hyb, ~] = g_all(t_hyb);
    is_f_hyb = all(gg_hyb <= 1e-4);
    
    % 3. Print the comparison
    fprintf('  Standard GSS: t* = %.4f mm, W* = %.6f kg, iterations: %d\n', ...
            t_gss, obj1D(t_gss), n_gss);
    fprintf('  Hybrid Search:t* = %.4f mm, W* = %.6f kg, feasible: %s, iterations: %d\n', ...
            t_hyb, W_hyb, mat2str(is_f_hyb), n_hyb);
            
    if exist('feas_rec', 'var') && any(feas_rec)
        fprintf('  fminsearch:   t* = %.4f mm, W* = %.6f kg\n', t_best_1d, bestW_1d);
        fprintf('  Difference |W_hyb - W_fms| = %.2e kg\n', abs(W_hyb - bestW_1d));
    end

    % -----------------------------------------------------------------------
    % Figure 10: Custom Algorithm Convergence
    % -----------------------------------------------------------------------
    figure('Name', 'Fig 10: Hybrid Algorithm Convergence', 'Position', [150 350 700 420]);
    
    semilogy(1:length(hist_hyb), hist_hyb, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('Iteration');
    ylabel('Bracket/Spread width (mm)');
    title('Custom Hybrid Search Convergence');
    grid on;
    
    hold on;
    switch_iter = find(hist_hyb <= 0.1 * hist_hyb(1), 1, 'first'); 
    if ~isempty(switch_iter)
        xline(switch_iter, 'r--', 'Switch to Parabolic', 'LabelVerticalAlignment', 'bottom');
    end
    hold off;

    % -----------------------------------------------------------------------
    % Figure 11: Side-by-Side Algorithm Comparison
    % -----------------------------------------------------------------------
    figure('Name', 'Fig 11: Algorithm Comparison', 'Position', [200 400 700 420]);
    
    semilogy(1:length(hist_gss), hist_gss, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4,'DisplayName', 'Standard GSS');
    hold on;
    semilogy(1:length(hist_hyb), hist_hyb, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Hybrid (GSS + Parabolic)');
    
    if ~isempty(switch_iter)
        xline(switch_iter, 'r:', 'Parabolic Stage Engages', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
    
    xlabel('Iteration Number', 'FontWeight','bold');
    ylabel('Interval Width (mm)', 'FontWeight','bold');
    title('Convergence Comparison: GSS vs. Hybrid Search');
    legend('Location', 'northeast');
    grid on;
    hold off;
    
    fprintf('  Section 4 complete.\n');
end