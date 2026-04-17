function section5_optimization(p)
% SECTION5_OPTIMIZATION  Full 4-D constrained optimization
%
% Method A: fmincon (SQP) 
% Method B: ALM + SD
% Method C: ALM + CG
% Method D: SLP with move limits


    lb = p.lb(:);
    ub = p.ub(:);

    obj = @(x) beam_weight(x, p);
    con = @(x) beam_constraints(x, p);

    %% Method A: Fmincon (Baseline)
    fprintf('\n  --- Method A: fmincon (SQP) ---\n');

    opts_fmincon = optimoptions('fmincon', ...
        'Algorithm',          'sqp', ...
        'Display',            'off', ...
        'OptimalityTolerance', 1e-10, ...
        'ConstraintTolerance', 1e-10, ...
        'StepTolerance',       1e-12, ...
        'MaxIterations',       2000);

    % Multiple starting points
    n_starts = 16;
    rng(42);
    x0_all = zeros(n_starts, 4);
    for i = 1:4
        x0_all(:, i) = lb(i) + rand(n_starts, 1) * (ub(i) - lb(i));
    end

    results_A = [];

   % Print Table Header
    fprintf('  %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s | %-10s | %-6s\n', ...
            'h0','b0','t0','th0', 'h*','b*','t*','th*', 'W* (kg)', 'Iters');
    fprintf('  %s\n', repmat('-', 1, 100));
    
    for i = 1:n_starts
        x0 = x0_all(i, :)';
        try
            [x_opt, W_opt, ef, output] = fmincon(obj, x0, [],[],[],[], lb, ub, ...
                @(x) deal(con(x), []), opts_fmincon);
                
            [c_opt, ~] = con(x_opt);
            
            % 1. Success: Converged AND Feasible
            if ef > 0 && all(c_opt <= 1e-4)
                results_A(end+1, :) = [x0', x_opt', W_opt, output.iterations];
                
                fprintf('  %-8.1f %-8.1f %-8.2f %-8.1f | %-8.2f %-8.2f %-8.3f %-8.2f | %-10.6f | %-6d\n', ...
                        x0(1), x0(2), x0(3), x0(4), x_opt(1), x_opt(2), x_opt(3), x_opt(4), W_opt, output.iterations);
                        
            % 2. Failure: Did not converge or resulted in an infeasible design
            else
                % Print the starting point on the left, and the failure reason on the right
                if ef <= 0
                    fprintf('  %-8.1f %-8.1f %-8.2f %-8.1f | [Run %2d] FAILED: Solver stopped (Exit Flag: %d)\n', ...
                            x0(1), x0(2), x0(3), x0(4), i, ef);
                else
                    fprintf('  %-8.1f %-8.1f %-8.2f %-8.1f | [Run %2d] FAILED: Converged to infeasible design\n', ...
                            x0(1), x0(2), x0(3), x0(4), i);
                end
            end
            
        catch
            % 3. Critical Failure: Code crashed
            fprintf('  %-8.1f %-8.1f %-8.2f %-8.1f | [Run %2d] FAILED: fmincon crashed mathematically\n', ...
                    x0(1), x0(2), x0(3), x0(4), i);
        end
    end

    if isempty(results_A)
        error('fmincon found no feasible solution. Check parameters.');
    end

    [W_opt_A, idx_A] = min(results_A(:, 9));
    x_opt_A = results_A(idx_A, 5:8)';
    iter_A = results_A(idx_A, 10);

    [~, ~, info_A] = beam_constraints(x_opt_A, p);
    fprintf('\n  fmincon optimum:\n');
    fprintf('    h*     = %.2f mm\n',   x_opt_A(1));
    fprintf('    b*     = %.2f mm\n',   x_opt_A(2));
    fprintf('    t*     = %.3f mm\n',   x_opt_A(3));
    fprintf('    theta* = %.2f deg\n',  x_opt_A(4));
    fprintf('    W*     = %.6f kg\n',   W_opt_A);
    fprintf('    sigma  = %.1f MPa  (limit %.0f)\n', info_A.sigma_bend, p.sigma_allow);
    fprintf('    tau    = %.1f MPa  (limit %.0f)\n', info_A.tau_torsion, p.tau_allow);
    fprintf('    delta  = %.2f mm   (limit %.0f)\n', info_A.delta, p.delta_max);
    fprintf('    phi    = %.2f deg  (limit %.1f)\n', info_A.phi_deg, p.phi_max*180/pi);

    %% Method C: Augmented Lagrangian
    
    % Small matrix of 3 carefully chosen starting points
    % Rows: [h, b, t, theta]
    x0_options = [
        65,  55, 2.0,  45;    % Point 1: The "Heavy/Safe" guess
        30.1, 40.1, 1.51, 45; % Point 2: The Mid way point
        60,  20, 1.0,  15     % Point 3: The "Tall/Narrow" guess
    ]';
    
    % --- VARIABLE NORMALIZATION SETUP (Lecture 7, Slide 42) ---
    % Translator 1: Normalized Space [0, 1] -> Physical Space [lb, ub]
    to_physical = @(xn) lb + xn .* (ub - lb);
    
    % Translator 2: Physical Space [lb, ub] -> Normalized Space [0, 1]
    to_normal = @(x) (x - lb) ./ (ub - lb);
    
    % Wrapped functions: Solver uses 'xn' (0 to 1), translated to physical limits
    obj_norm = @(xn) beam_weight(to_physical(xn), p);
    con_norm = @(xn) beam_constraints(to_physical(xn), p);
    
    % The inner solver bounds are now purely 0 and 1
    lb_norm = zeros(4, 1);
    ub_norm = ones(4, 1);
    
    fprintf('\n  --- Method C: ALM + SD & ALM + CG ---\n');
    
    % Loop through your safe starting points
    for i = 1:size(x0_options, 2)
        % Get physical start point and immediately normalize it
        x_start_phys = x0_options(:, i);
        x_start_norm = to_normal(x_start_phys); 
        
        fprintf('\n  Testing ALM from Start Point %d...\n', i);
        
        % Setup Algorithm Options & ALM Parameters
        max_al_iter   = 40;
        p_weight_init = 1.0; 
        
        % Options for inner solvers
        sd_options.TolF     = 1e-7; % Loosened for early ALM iterations 
        sd_options.MaxIter  = 1000;
        sd_options.h        = 1e-5;
        sd_options.beta     = 0.0;  % SD has no momentum
        
        cg_options          = sd_options;
        
        % -----------------------------------------------------------
        % Run ALM + Steepest Descent (SD)
        % -----------------------------------------------------------
        x_al_sd       = x_start_norm;
        
        % Evaluate constraints once to capture exact shape
        [g_init_sd, ~] = con_norm(x_al_sd);
        lambda_sd      = zeros(size(g_init_sd)); 
        
        p_weight_sd   = p_weight_init;
        iter_total_sd = 0;

        % Initialize the full history array before the loop
        full_hist_sd = [];
        
        for k = 1:max_al_iter
            % Objective in 0-to-1 space with current lambda
            obj_al_sd = @(xn) obj_norm(xn) + ...
                p_weight_sd * sum(max(0, con_norm(xn) + lambda_sd / (2 * p_weight_sd)).^2);
            
            % Run inner solver with normalized bounds
            [x_al_sd, ~, hist_sd] = steepest_descent(obj_al_sd, x_al_sd, lb_norm, ub_norm, sd_options);
            iter_total_sd = iter_total_sd + size(hist_sd, 1) - 1;

            % Append the iteration's history to the master record
            full_hist_sd = [full_hist_sd; hist_sd];
            
            % Evaluate constraints
            [g_sd, ~] = con_norm(x_al_sd);
            viol_sd = max(0, g_sd);
            max_viol_sd = max(viol_sd);
            
            % NEW: Conditional Penalty Update
            if k == 1 || max_viol_sd > 0.25 * max_viol_prev_sd
                p_weight_sd = p_weight_sd * 2.0; % Only increase if not improving fast enough
            end
            max_viol_prev_sd = max_viol_sd;
            
            % Update multipliers
            lambda_sd = lambda_sd + 2 * p_weight_sd * viol_sd;
            
        end
        
        % -----------------------------------------------------------
        % Run ALM + Conjugate Gradient (CG)
        % -----------------------------------------------------------
        x_al_cg       = x_start_norm;
        
        % FIX: Evaluate constraints once to capture exact shape
        [g_init_cg, ~] = con_norm(x_al_cg);
        lambda_cg      = zeros(size(g_init_cg));
        
        p_weight_cg   = p_weight_init;
        iter_total_cg = 0;

        % Initialize the full history array before the loop
        full_hist_cg = [];
        
        for k = 1:max_al_iter
            % Objective in 0-to-1 space with current lambda
            obj_al_cg = @(xn) obj_norm(xn) + ...
                p_weight_cg * sum(max(0, con_norm(xn) + lambda_cg / (2 * p_weight_cg)).^2);
            
            % Run inner solver with normalized bounds
            [x_al_cg, ~, hist_cg] = conjugate_gradient(obj_al_cg, x_al_cg, lb_norm, ub_norm, cg_options);
            iter_total_cg = iter_total_cg + size(hist_cg, 1) - 1;

            % Append the iteration's history to the master record
            full_hist_cg = [full_hist_cg; hist_cg];

            % Evaluate constraints
            [g_cg, ~] = con_norm(x_al_cg);
            viol_cg = max(0, g_cg);
            max_viol_cg = max(viol_cg);
            
            % Conditional Penalty Update
            if k == 1 || max_viol_cg > 0.25 * max_viol_prev_cg
                p_weight_cg = p_weight_cg * 2.0; % Only increase if not improving fast enough
            end
            max_viol_prev_cg = max_viol_cg;
            
            % Update multipliers
            lambda_cg = lambda_cg + 2 * p_weight_cg * viol_cg;
        end
        
        % -----------------------------------------------------------
        % Compare Results
        % -----------------------------------------------------------
        % Translate the final 0-to-1 optimal points back to physical reality
        x_al_sd_final = to_physical(x_al_sd);
        x_al_cg_final = to_physical(x_al_cg);
        
        % Calculate physical weights
        W_sd = beam_weight(x_al_sd_final, p);
        W_cg = beam_weight(x_al_cg_final, p);
        
        fprintf('    ALM + SD: h*=%.2f  b*=%.2f  t*=%.3f  theta*=%.2f | W* = %.6f kg (%d inner iter)\n', ...
                x_al_sd_final(1), x_al_sd_final(2), x_al_sd_final(3), x_al_sd_final(4), W_sd, iter_total_sd);
                
        fprintf('    ALM + CG: h*=%.2f  b*=%.2f  t*=%.3f  theta*=%.2f | W* = %.6f kg (%d inner iter)\n', ...
                x_al_cg_final(1), x_al_cg_final(2), x_al_cg_final(3), x_al_cg_final(4), W_cg, iter_total_cg);
    end
    
    % Store the final physical points for the comparison table below
    x_al_sd = x_al_sd_final;
    x_al_cg = x_al_cg_final;

    %% Method D: SLP with Move Limits (Multi-Start)
    fprintf('\n  --- Method D: SLP with Move Limits (Multi-Start) ---\n');
    ml = [5; 4; 0.3; 10];   % move limits for [h, b, t, theta]
    slp_options.TolX    = 1e-8;
    slp_options.MaxIter = 200;
    slp_options.h       = 1e-5;
    
    % 1. Setup Multiple Starting Points (Using the same logic as fmincon)
    n_starts_slp = 16;
    rng(42); % Use the same seed for a fair comparison
    x0_slp_all = zeros(n_starts_slp, 4);
    for i = 1:4
        x0_slp_all(:, i) = lb(i) + rand(n_starts_slp, 1) * (ub(i) - lb(i));
    end
    
    results_D = [];
    best_hist_slp = []; % To keep track of the path for the absolute best run
    
    fprintf('  %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s | %-10s | %-6s\n', ...
            'h0','b0','t0','th0', 'h*','b*','t*','th*', 'W* (kg)', 'Iters');
    fprintf('  %s\n', repmat('-', 1, 100));
    
    % 2. Run SLP from every starting point
    for i = 1:n_starts_slp
        x0_current = x0_slp_all(i, :)';
        try
            [x_slp, W_slp, current_hist_slp] = slp_move_limits(obj, ...
                @(x) deal(con(x), []), x0_current, lb, ub, ml, slp_options);
                
            [c_slp, ~] = con(x_slp);
            
            % 3. Check for Feasibility
            if all(c_slp <= 1e-4)
                iter_current = size(current_hist_slp, 1) - 1;
                results_D(end+1, :) = [x0_current', x_slp', W_slp, iter_current];
                
                fprintf('  %-8.1f %-8.1f %-8.2f %-8.1f | %-8.2f %-8.2f %-8.3f %-8.2f | %-10.6f | %-6d\n', ...
                        x0_current(1), x0_current(2), x0_current(3), x0_current(4), ...
                        x_slp(1), x_slp(2), x_slp(3), x_slp(4), W_slp, iter_current);
                        
                if W_slp == min(results_D(:, 9))
                    best_hist_slp = current_hist_slp;
                end
            else
                % NEW: Print if it converged but violated constraints
                fprintf('  [Run %2d] FAILED: Converged to an infeasible design.\n', i);
            end
        catch
            % NEW: Print if linprog completely crashed
            fprintf('  [Run %2d] FAILED: linprog crashed (LP Infeasible).\n', i);
        end
    end
    
    if isempty(results_D)
        error('SLP found no feasible solution from any starting point.');
    end
    
    % 4. Extract the absolute best overall run
    [W_opt_D, idx_D] = min(results_D(:, 9));
    x_opt_D = results_D(idx_D, 5:8)';
    iter_SLP = results_D(idx_D, 10);
    hist_slp = best_hist_slp; % Restore the best history so your later plots work perfectly
    
    [c_slp_final, ~] = con(x_opt_D);
    
    fprintf('\n  Best SLP optimum found:\n');
    fprintf('    h* = %.2f, b* = %.2f, t* = %.3f, theta* = %.2f\n', x_opt_D(1), x_opt_D(2), x_opt_D(3), x_opt_D(4));
    fprintf('    W* = %.6f kg\n', W_opt_D);
    fprintf('    Iterations: %d\n', iter_SLP);
    fprintf('    Constraints: [%.4f, %.4f, %.4f, %.4f, %.4f]\n', c_slp_final);

    % =======================================================================
    % FINAL CONSOLE REPORT: Comparison Table
    % =======================================================================
    fprintf('\n  ===== Comparison of All Methods =====\n');
    fprintf('  %-20s %-7s %-7s %-8s %-9s %-10s %-6s\n', ...
            'Method', 'h*', 'b*', 't*', 'theta*', 'W* (kg)', 'Iters');
    fprintf('  %s\n', repmat('-', 1, 75));
    
    % Method A: fmincon
    fprintf('  %-20s %-7.2f %-7.2f %-8.5f %-9.4f %-10.6f %-6d\n', ...
            'fmincon (SQP)', x_opt_A(1), x_opt_A(2), x_opt_A(3), x_opt_A(4), W_opt_A, iter_A);
            
    % Method B: ALM + Steepest Descent
    fprintf('  %-20s %-7.2f %-7.2f %-8.3f %-9.2f %-10.6f %-6d\n', ...
            'ALM + SD', x_al_sd(1), x_al_sd(2), x_al_sd(3), x_al_sd(4), W_sd, iter_total_sd);
            
    % Method C: ALM + Conjugate Gradient
    fprintf('  %-20s %-7.2f %-7.2f %-8.3f %-9.2f %-10.6f %-6d\n', ...
            'ALM + CG', x_al_cg(1), x_al_cg(2), x_al_cg(3), x_al_cg(4), W_cg, iter_total_cg);
            
    % Method D: SLP
    fprintf('  %-20s %-7.2f %-7.2f %-8.3f %-9.2f %-10.6f %-6d\n', ...
            'SLP move limits', x_opt_D(1), x_opt_D(2), x_opt_D(3), x_opt_D(4), W_opt_D, iter_SLP);
    fprintf('\n');
    
    % -----------------------------------------------------------------------
    % 3D PLOT: ALM Landscape Warping (Effect of Initial p_weight)
    % -----------------------------------------------------------------------
    % fprintf('\n  Generating Figure 16: ALM Landscape Warping...\n');
    % 
    % % 1. Setup the 2D Slice Grid
    % N_grid = 125;
    % h_vec = linspace(p.lb(1), p.ub(1), N_grid);
    % b_vec = linspace(p.lb(2), p.ub(2), N_grid);
    % [H_grid, B_grid] = meshgrid(h_vec, b_vec);
    % 
    % % Fix the other two variables
    % t_fix = 0.547; 
    % th_fix = 22;
    % 
    % W_grid = zeros(N_grid, N_grid);
    % P_grid = zeros(N_grid, N_grid); % The Penalty Term
    % G_grid = zeros(N_grid, N_grid, 5); % Constraint term
    % 
    % % 2. Pre-calculate the Physics over the Grid
    % for ii = 1:N_grid
    %     for jj = 1:N_grid
    %         x_current = [H_grid(ii,jj); B_grid(ii,jj); t_fix; th_fix];
    % 
    %         % Get Objective (Weight)
    %         W_grid(ii,jj) = beam_weight(x_current, p);
    % 
    %         % Get Constraints
    %         [g_vals, ~] = beam_constraints(x_current, p);
    % 
    %         % Save individual constraints for plotting later
    %         G_grid(ii,jj,:) = g_vals; % <--- NEW: Save the constraints
    % 
    %         % Calculate the squared penalty sum
    %         P_grid(ii,jj) = sum(max(0, g_vals).^2);
    %     end
    % end
    % 
    % % 3. Create the Visualization (Updated for 3 weights)
    % % Made the figure wider and shorter to fit a 1x3 side-by-side layout perfectly
    % figure('Name', 'ALM Landscape Warping', 'Position', [100, 200, 1500, 500]);
    % 
    % % NEW: Only testing 3 weights
    % test_weights = [1, 50, 500];
    % 
    % % NEW: Loop exactly 3 times
    % for i = 1:3
    %     pw = test_weights(i);
    % 
    %     % The ALM Math: F = Objective + p_weight * Penalty
    %     F_Total = W_grid + pw * P_grid;
    % 
    %     % Find the absolute lowest point in this specific landscape
    %     [min_val, min_idx] = min(F_Total(:));
    %     [row, col] = ind2sub(size(F_Total), min_idx);
    %     min_h = H_grid(row, col);
    %     min_b = B_grid(row, col);
    %     min_weight = W_grid(row, col);
    % 
    %     % NEW: Use a 1x3 subplot grid instead of 2x2
    %     subplot(1, 3, i);
    %     surf(H_grid, B_grid, F_Total, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    %     hold on;
    % 
    %     % Add contour lines on the floor to show the "valley" shape
    %     contour3(H_grid, B_grid, F_Total, 30, 'k');
    % 
    %     % Plot the temporary "minimum" the solver will be attracted to
    %     plot3(min_h, min_b, min_val, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    % 
    %     % ---------------------------------------------------------
    %     % DRAW THE EXACT PHYSICAL CONSTRAINTS ON THE FLOOR
    %     % ---------------------------------------------------------
    %     % g1: Bending Stress (Red line)
    %     contour3(H_grid, B_grid, G_grid(:,:,1), [0 0], 'r', 'LineWidth', 2.5);
    % 
    %     % g3: Tip Deflection (Black line)
    %     contour3(H_grid, B_grid, G_grid(:,:,3), [0 0], 'k', 'LineWidth', 2.5);
    % 
    %     % g4: Twist Angle (Blue line)
    %     contour3(H_grid, B_grid, G_grid(:,:,4), [0 0], 'b', 'LineWidth', 2.5);
    % 
    %     % Aesthetics
    %     colormap(parula);
    %     view(50, 20); 
    %     xlabel('h (mm)', 'FontWeight', 'bold');
    %     ylabel('b (mm)', 'FontWeight', 'bold');
    %     zlabel('Total ALM Objective F(x)', 'FontWeight', 'bold');
    %     title({sprintf('Initial p\\_weight = %d', pw), ...
    %            sprintf('Temp Optimum: W = %.3f kg (Total F = %.2f)', min_weight, min_val)}, ...
    %            'FontSize', 11, 'FontWeight', 'bold');
    % 
    %     % Lock the Z-axis so the warping effect is visually obvious across all subplots
    %     zlim([0, max(F_Total(:))]); 
    %     grid on; hold off;
    % end
    % 
    % % =======================================================================
    % % REPORT PLOT: fmincon vs. ALM + SD vs. ALM + CG
    % % =======================================================================
    % fprintf('\n  Generating Plot: SD vs CG Optimization Paths...\n');
    % 
    % % 1. Setup the 2D Slice Grid
    % % We fix thickness and ply angle to their final optimal values
    % t_fix = 0.54748;  
    % th_fix = 21.31917; 
    % 
    % N_grid = 100;
    % h_vec = linspace(p.lb(1), p.ub(1), N_grid);
    % b_vec = linspace(p.lb(2), p.ub(2), N_grid);
    % [H_grid, B_grid] = meshgrid(h_vec, b_vec);
    % 
    % W_grid = zeros(N_grid, N_grid);
    % G_grid = zeros(N_grid, N_grid, 5); 
    % 
    % % 2. Pre-calculate Weight and Constraints over the Grid
    % for ii = 1:N_grid
    %     for jj = 1:N_grid
    %         x_eval = [H_grid(ii,jj); B_grid(ii,jj); t_fix; th_fix];
    %         W_grid(ii,jj) = beam_weight(x_eval, p);
    %         [g_vals, ~] = beam_constraints(x_eval, p);
    %         G_grid(ii,jj,:) = g_vals; 
    %     end
    % end
    % 
    % % 3. Create the Visualization
    % figure('Name', 'SD vs CG Trajectories', 'Position', [150, 150, 800, 600]);
    % hold on;
    % 
    % % --- BACKGROUND: Weight Contours ---
    % contour(H_grid, B_grid, W_grid, 25, 'DisplayName', 'W (kg)');
    % colormap("parula")
    % 
    % % --- CONSTRAINTS: True Physical Boundaries (g_i = 0) ---
    % % 1. Bending Stress (Red)
    % contour(H_grid, B_grid, G_grid(:,:,1), [0 0], 'r', 'LineWidth', 2.5, 'DisplayName', 'Bending Limit');
    % 
    % % 2. Shear Stress (Green)
    % contour(H_grid, B_grid, G_grid(:,:,2), [0 0], 'Color', [0 0.8 0], 'LineWidth', 2.5, 'DisplayName', 'Shear Limit', 'HandleVisibility', 'off');
    % 
    % % 3. Deflection (Black)
    % contour(H_grid, B_grid, G_grid(:,:,3), [0 0], 'k', 'LineWidth', 2.5, 'DisplayName', 'Deflection Limit', 'HandleVisibility', 'off');
    % 
    % % 4. Twist Angle (Blue)
    % contour(H_grid, B_grid, G_grid(:,:,4), [0 0], 'b', 'LineWidth', 2.5, 'DisplayName', 'Twist Limit');
    % 
    % % 5. Buckling (Teal)
    % contour(H_grid, B_grid, G_grid(:,:,5), [0 0], 'Color', [0 1 1], 'LineWidth', 2.5, 'DisplayName', 'Buckling Limit');
    % 
    % 
    % % -------------------------------------------------------------------------
    % % Convert Normalized Histories back to Physical Space for plotting
    % % -------------------------------------------------------------------------
    % phys_hist_sd = full_hist_sd;
    % phys_hist_cg = full_hist_cg;
    % 
    % % Scale ALL 4 variables back to physical space using your lb and ub arrays
    % phys_hist_sd(:, 2) = lb(1) + full_hist_sd(:, 2) .* (ub(1) - lb(1)); % h
    % phys_hist_sd(:, 3) = lb(2) + full_hist_sd(:, 3) .* (ub(2) - lb(2)); % b
    % phys_hist_sd(:, 4) = lb(3) + full_hist_sd(:, 4) .* (ub(3) - lb(3)); % t
    % phys_hist_sd(:, 5) = lb(4) + full_hist_sd(:, 5) .* (ub(4) - lb(4)); % theta
    % 
    % phys_hist_cg(:, 2) = lb(1) + full_hist_cg(:, 2) .* (ub(1) - lb(1)); % h
    % phys_hist_cg(:, 3) = lb(2) + full_hist_cg(:, 3) .* (ub(2) - lb(2)); % b
    % phys_hist_cg(:, 4) = lb(3) + full_hist_cg(:, 4) .* (ub(3) - lb(3)); % t
    % phys_hist_cg(:, 5) = lb(4) + full_hist_cg(:, 5) .* (ub(4) - lb(4)); % theta
    % 
    % % -------------------------------------------------------------------------
    % % 3a. Capture fmincon (SQP) path from the shared starting point
    % % -------------------------------------------------------------------------
    % global hist_fmincon;
    % hist_fmincon = []; % Initialize empty history array
    % 
    % % Reverted to your preferred method: Extract the exact physical start point
    % x0_shared = [phys_hist_sd(1,2); phys_hist_sd(1,3); phys_hist_sd(1,4); phys_hist_sd(1,5)];
    % 
    % % Setup fmincon to use a custom Output Function to track its path
    % opts_trace = optimoptions('fmincon', ...
    %     'Algorithm', 'sqp', ...
    %     'Display', 'off', ...
    %     'OutputFcn', @outfun_fmincon);
    % 
    % % Run fmincon once just to capture the history
    % [x_fmincon_final, ~] = fmincon(obj, x0_shared, [],[],[],[], lb, ub, ...
    %                                @(x) deal(con(x), []), opts_trace);
    % 
    % % Scale h (Column 2) back to physical
    % phys_hist_sd(:, 2) = p.lb(1) + full_hist_sd(:, 2) .* (p.ub(1) - p.lb(1));
    % phys_hist_cg(:, 2) = p.lb(1) + full_hist_cg(:, 2) .* (p.ub(1) - p.lb(1));
    % 
    % % Scale b (Column 3) back to physical
    % phys_hist_sd(:, 3) = p.lb(2) + full_hist_sd(:, 3) .* (p.ub(2) - p.lb(2));
    % phys_hist_cg(:, 3) = p.lb(2) + full_hist_cg(:, 3) .* (p.ub(2) - p.lb(2));
    % 
    % % --- THE PATHS: Plot using the new physical coordinates ---
    % % Plot Steepest Descent (SD)
    % plot(phys_hist_sd(:,2), phys_hist_sd(:,3), 'k.-', 'LineWidth', 1.0, 'MarkerSize', 8, ...
    %      'DisplayName', 'ALM + SD (Zig-Zag)');
    % 
    % % Plot Conjugate Gradient (CG)
    % plot(phys_hist_cg(:,2), phys_hist_cg(:,3), 'g.--', 'LineWidth', 1.0, 'MarkerSize', 10, ...
    %      'DisplayName', 'ALM + CG (Momentum)');
    % 
    % % Plot fmincon (SQP)
    % plot(hist_fmincon(:,2), hist_fmincon(:,3), 'm.--', 'LineWidth', 1.0, 'MarkerSize', 10, ...
    %      'DisplayName', 'fmincon (SQP)');
    % 
    % % Mark the Start and End Points
    % plot(phys_hist_sd(1,2), phys_hist_sd(1,3), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Start Point');
    % plot(x_al_sd_final(1), x_al_sd_final(2), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'DisplayName', 'Optimum');
    % 
    % % Formatting
    % xlabel('Height h (mm)', 'FontSize', 12, 'FontWeight', 'bold'); 
    % ylabel('Width b (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % title('fmincon vs. ALM + SD vs. ALM + CG in h-b design space', 'FontSize', 14, 'FontWeight', 'bold');
    % legend('Location', 'eastoutside', 'FontSize', 11);
    % 
    % % Dynamically set limits based on the physical paths so the zig-zag is visible
    % xlim([min([phys_hist_sd(:,2); phys_hist_cg(:,2)]) - 2, max([phys_hist_sd(:,2); phys_hist_cg(:,2)]) + 2]);
    % ylim([min([phys_hist_sd(:,3); phys_hist_cg(:,3)]) - 2, max([phys_hist_sd(:,3); phys_hist_cg(:,3)]) + 2]);
    % 
    % grid on; hold off;
    % 
    % % =======================================================================
    % % REPORT PLOT: Global SLP Linearization of ALL Constraints (Clean)
    % % =======================================================================
    % fprintf('\n  Generating Plot: Clean Linearization of All Constraints...\n');
    % 
    % % 1. Set the specific design space slice requested
    % t_plot  = 0.547;  
    % th_plot = 21.32;  
    % 
    % % Mathematical evaluation point for the Taylor Series (not plotted)
    % hk = 70; 
    % bk = 60; 
    % xk = [hk; bk; t_plot; th_plot];
    % 
    % % 2. Create the Grid
    % N_grid = 150;
    % h_vec = linspace(25, 75, N_grid); 
    % b_vec = linspace(15, 65, N_grid); 
    % [H_grid, B_grid] = meshgrid(h_vec, b_vec);
    % 
    % % Pre-allocate constraint matrices
    % G_true = zeros(N_grid, N_grid, 5);
    % G_lin  = zeros(N_grid, N_grid, 5);
    % W_grid = zeros(N_grid, N_grid);
    % 
    % % 3. Calculate True Constraints & Weight over the grid
    % for ii = 1:N_grid
    %     for jj = 1:N_grid
    %         x_eval = [H_grid(ii,jj); B_grid(ii,jj); t_plot; th_plot];
    % 
    %         % NEW: Calculate Objective Function
    %         W_grid(ii,jj) = beam_weight(x_eval, p);
    % 
    %         % Existing Constraint Calculation
    %         [g_vals, ~] = beam_constraints(x_eval, p);
    %         G_true(ii,jj,:) = g_vals; 
    %     end
    % end
    % 
    % % 4. Calculate Gradients at the Evaluation Point (xk)
    % h_step = slp_options.h;
    % [gk_vals, ~] = beam_constraints(xk, p); 
    % 
    % xk_h = xk; xk_h(1) = xk_h(1) + h_step;
    % [g_h, ~] = beam_constraints(xk_h, p);
    % dg_dh = (g_h - gk_vals) / h_step; 
    % 
    % xk_b = xk; xk_b(2) = xk_b(2) + h_step;
    % [g_b, ~] = beam_constraints(xk_b, p);
    % dg_db = (g_b - gk_vals) / h_step; 
    % 
    % % 5. Calculate Linearized Constraints over the grid
    % for c = 1:5
    %     G_lin(:,:,c) = gk_vals(c) + dg_dh(c)*(H_grid - hk) + dg_db(c)*(B_grid - bk);
    % end
    % 
    % % 6. Create the Visualization
    % figure('Name', 'Clean SLP Linearization', 'Position', [250, 150, 900, 700]);
    % hold on;
    % 
    % % NEW: Plot the Linear Objective Function Contours in light grey
    % contour(H_grid, B_grid, W_grid, 20, ...
    %         'ShowText', 'on', 'LabelSpacing', 400, 'DisplayName', 'W (kg)');
    % colormap('parula');
    % 
    % % Custom Color Mapping [R, G, B]
    % custom_colors = [
    %     1, 0, 0;        % 1: Bending Stress (Red)
    %     0, 1, 0;      % 2: Shear Stress (Green - slightly darker for visibility)
    %     0, 0, 0;        % 3: Deflection (Black)
    %     0, 0, 1;        % 4: Twist Angle (Blue)
    %     0, 1, 1     % 5: Buckling (Teal)
    % ];
    % 
    % con_names = {'Bending Stress', 'Shear Stress', 'Deflection', 'Twist Angle', 'Buckling'};
    % 
    % % Plot ALL True and Linearized Constraints
    % for c = 1:5
    %     % True Constraint (Solid Line, thick)
    %     contour(H_grid, B_grid, G_true(:,:,c), [0 0], 'Color', custom_colors(c,:), ...
    %             'LineWidth', 2.5, 'DisplayName', [con_names{c}]);
    % 
    %     % Linearized Constraint (Dashed Line, matching color)
    %     contour(H_grid, B_grid, G_lin(:,:,c), [0 0], '--', 'Color', custom_colors(c,:), ...
    %             'LineWidth', 2.0, 'HandleVisibility', 'off'); 
    % end
    % 
    % % Plot the optimum
    % plot(hk, bk, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'DisplayName', 'Optimum');
    % 
    % % Formatting
    % xlabel('Height h (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('Width b (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % title(sprintf('True vs Linearized Constraints (t = %.2f mm, \\theta = %.2f^\\circ)', t_plot, th_plot), ...
    %       'FontSize', 14, 'FontWeight', 'bold');
    % 
    % legend('Location', 'eastoutside', 'FontSize', 11);
    % xlim([25, 75]);
    % ylim([15, 65]);
    % grid on; hold off;
    % 
    % % -----------------------------------------------------------------------
    % % Figure: Sensitivity analysis — W* vs tip load P
    % % -----------------------------------------------------------------------
    % fprintf('\n  --- Sensitivity Analysis: W* vs Tip Load P ---\n');
    % 
    % P_list = linspace(500, 4000, 20);
    % W_sens = zeros(size(P_list));
    % x_sens = zeros(length(P_list), 4);
    % 
    % for i = 1:length(P_list)
    %     p_temp   = p;
    %     p_temp.P = P_list(i);
    %     obj_t = @(x) beam_weight(x, p_temp);
    %     con_t = @(x) beam_constraints(x, p_temp);
    %     x0    = x_opt_A;
    %     try
    %         [x_s, W_s, ef] = fmincon(obj_t, x0, [],[],[],[], lb, ub, ...
    %             @(x) deal(con_t(x), []), opts_fmincon);
    %         [c_s, ~] = con_t(x_s);
    %         if ef > 0 && all(c_s <= 1e-4)
    %             W_sens(i)    = W_s;
    %             x_sens(i, :) = x_s';
    %         else
    %             W_sens(i) = NaN;
    %         end
    %     catch
    %         W_sens(i) = NaN;
    %     end
    % end
    % 
    % figure('Name', 'Sensitivity Analysis', 'Position', [150 80 1000 400]);
    % 
    % subplot(1,2,1);
    % plot(P_list, W_sens*1000, 'b-', 'LineWidth', 1.5);
    % hold on;
    % xline(p.P, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Nominal P');
    % xlabel('Tip load P (N)', 'FontWeight', 'bold');
    % ylabel('Optimal weight W* (g)', 'FontWeight', 'bold');
    % title('Optimal weight vs tip load');
    % legend({'W*(P)', 'Assumed Tip Load'}); grid on; hold off;
    % 
    % subplot(1,2,2);
    % plot(P_list, x_sens(:,1), 'r-',  'LineWidth', 1.5, 'DisplayName', 'h* (mm)'); hold on;
    % plot(P_list, x_sens(:,2), 'g-',  'LineWidth', 1.5, 'DisplayName', 'b* (mm)');
    % plot(P_list, x_sens(:,3)*10, 'k-', 'LineWidth', 1.5, 'DisplayName', 't* x10 (mm)');
    % plot(P_list, x_sens(:,4), 'b-',  'LineWidth', 1.5, 'DisplayName', '\theta* (deg)');
    % xline(p.P, 'k--', 'LineWidth', 1, 'DisplayName', 'Assumed Tip Load');
    % xlabel('Tip load P (N)', 'FontWeight', 'bold');
    % ylabel('Optimal variable value', 'FontWeight', 'bold');
    % title('Optimal variables vs load');
    % legend(); grid on; hold off;
    % 
    % % =======================================================================
    % % REPORT PLOT: 3D Design Spaces (h-b, h-t, b-t) with Constraints
    % % =======================================================================
    % fprintf('\n  Generating Plot: 3D Design Space Slices with Constraints...\n');
    % 
    % % 1. Setup Fixed Parameters (using the optimums)
    % opt_h  = 70;
    % opt_b  = 60;
    % opt_t  = 0.547;
    % opt_th = 21.32;
    % opt_W  = 0.204976;
    % 
    % N_grid = 50;
    % h_vec = linspace(25, 75, N_grid);
    % b_vec = linspace(15, 65, N_grid);
    % t_vec = linspace(0.5, 3.0, N_grid);
    % 
    % % 2. Setup Figure (Made wider to fit 3 subplots perfectly)
    % figure('Name', '3D Weight Surfaces', 'Position', [50, 200, 1800, 600]);
    % 
    % % =========================================================
    % % SUBPLOT 1: Height (h) vs Width (b)
    % % =========================================================
    % subplot(1, 3, 1);
    % [H_grid1, B_grid1] = meshgrid(h_vec, b_vec);
    % W_grid_hb = zeros(N_grid, N_grid);
    % G_grid_hb = zeros(N_grid, N_grid, 5);
    % 
    % for ii = 1:N_grid
    %     for jj = 1:N_grid
    %         % Evaluate [h, b, t_fixed, theta_fixed]
    %         x_eval = [H_grid1(ii,jj); B_grid1(ii,jj); opt_t; opt_th];
    %         W_grid_hb(ii,jj) = beam_weight(x_eval, p);
    %         [g_vals, ~] = beam_constraints(x_eval, p);
    %         G_grid_hb(ii,jj,:) = g_vals;
    %     end
    % end
    % 
    % surf(H_grid1, B_grid1, W_grid_hb, 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    % hold on;
    % contour3(H_grid1, B_grid1, W_grid_hb, 20, 'k', 'HandleVisibility', 'off');
    % 
    % % Floor Constraints
    % contour3(H_grid1, B_grid1, G_grid_hb(:,:,1), [0 0], 'r', 'LineWidth', 2.5, 'DisplayName', 'Bending Limit');
    % contour3(H_grid1, B_grid1, G_grid_hb(:,:,2), [0 0], 'g', 'LineWidth', 2.5, 'DisplayName', 'Shear Limit');
    % contour3(H_grid1, B_grid1, G_grid_hb(:,:,3), [0 0], 'k', 'LineWidth', 2.5, 'DisplayName', 'Deflection Limit');
    % contour3(H_grid1, B_grid1, G_grid_hb(:,:,4), [0 0], 'b', 'LineWidth', 2.5, 'DisplayName', 'Twist Limit');
    % contour3(H_grid1, B_grid1, G_grid_hb(:,:,5), [0 0], 'Color', [0 1 1], 'LineWidth', 2.5, 'DisplayName', 'Buckling Limit');
    % 
    % % Optimum Star and Drop Line
    % plot3(opt_h, opt_b, opt_W, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', 'DisplayName', 'Optimum');
    % plot3([opt_h, opt_h], [opt_b, opt_b], [0, opt_W], 'r:', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % 
    % view(45, 20);
    % zlim([0, max(W_grid_hb(:))]);
    % xlabel('h (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('b (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % zlabel('W (kg)', 'FontSize', 12, 'FontWeight', 'bold');
    % title(sprintf('h vs b Design Space\n(Fixed t=%.3f, \\theta=%.2f^\\circ)', opt_t, opt_th), 'FontSize', 14, 'FontWeight', 'bold');
    % grid on; hold off;
    % 
    % % =========================================================
    % % SUBPLOT 2: Height (h) vs Thickness (t)
    % % =========================================================
    % subplot(1, 3, 2);
    % [H_grid2, T_grid1] = meshgrid(h_vec, t_vec);
    % W_grid_ht = zeros(N_grid, N_grid);
    % G_grid_ht = zeros(N_grid, N_grid, 5);
    % 
    % for ii = 1:N_grid
    %     for jj = 1:N_grid
    %         % Evaluate [h, b_fixed, t, theta_fixed]
    %         x_eval = [H_grid2(ii,jj); opt_b; T_grid1(ii,jj); opt_th];
    %         W_grid_ht(ii,jj) = beam_weight(x_eval, p);
    %         [g_vals, ~] = beam_constraints(x_eval, p);
    %         G_grid_ht(ii,jj,:) = g_vals;
    %     end
    % end
    % 
    % surf(H_grid2, T_grid1, W_grid_ht, 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    % hold on;
    % contour3(H_grid2, T_grid1, W_grid_ht, 20, 'k', 'HandleVisibility', 'off');
    % 
    % % Floor Constraints
    % contour3(H_grid2, T_grid1, G_grid_ht(:,:,1), [0 0], 'r', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % contour3(H_grid2, T_grid1, G_grid_ht(:,:,2), [0 0], 'Color', [0 0.8 0], 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % contour3(H_grid2, T_grid1, G_grid_ht(:,:,3), [0 0], 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % contour3(H_grid2, T_grid1, G_grid_ht(:,:,4), [0 0], 'b', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % contour3(H_grid2, T_grid1, G_grid_ht(:,:,5), [0 0], 'Color', [0 1 1], 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % 
    % plot3(opt_h, opt_t, opt_W, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', 'HandleVisibility', 'off');
    % plot3([opt_h, opt_h], [opt_t, opt_t], [0, opt_W], 'r:', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % 
    % colormap(parula);
    % view(45, 20);
    % zlim([0, max(W_grid_ht(:))]);
    % xlabel('h (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('t (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % zlabel('W (kg)', 'FontSize', 12, 'FontWeight', 'bold');
    % title(sprintf('h vs t Design Space\n(Fixed b=%.0f, \\theta=%.2f^\\circ)', opt_b, opt_th), 'FontSize', 14, 'FontWeight', 'bold');
    % 
    % grid on; hold off;
    % 
    % % =========================================================
    % % SUBPLOT 3: Width (b) vs Thickness (t)
    % % =========================================================
    % subplot(1, 3, 3);
    % [B_grid2, T_grid2] = meshgrid(b_vec, t_vec);
    % W_grid_bt = zeros(N_grid, N_grid);
    % G_grid_bt = zeros(N_grid, N_grid, 5);
    % 
    % for ii = 1:N_grid
    %     for jj = 1:N_grid
    %         % Evaluate [h_fixed, b, t, theta_fixed]
    %         x_eval = [opt_h; B_grid2(ii,jj); T_grid2(ii,jj); opt_th];
    %         W_grid_bt(ii,jj) = beam_weight(x_eval, p);
    %         [g_vals, ~] = beam_constraints(x_eval, p);
    %         G_grid_bt(ii,jj,:) = g_vals;
    %     end
    % end
    % 
    % surf(B_grid2, T_grid2, W_grid_bt, 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
    % hold on;
    % contour3(B_grid2, T_grid2, W_grid_bt, 20, 'k', 'HandleVisibility', 'off');
    % 
    % % Floor Constraints
    % contour3(B_grid2, T_grid2, G_grid_bt(:,:,1), [0 0], 'r', 'LineWidth', 2.5, 'DisplayName', 'Bending Limit');
    % contour3(B_grid2, T_grid2, G_grid_bt(:,:,2), [0 0], 'g' , 'LineWidth', 2.5, 'DisplayName', 'Shear Limit');
    % contour3(B_grid2, T_grid2, G_grid_bt(:,:,3), [0 0], 'k', 'LineWidth', 2.5, 'DisplayName', 'Deflection Limit');
    % contour3(B_grid2, T_grid2, G_grid_bt(:,:,4), [0 0], 'b', 'LineWidth', 2.5, 'DisplayName', 'Twist Limit');
    % contour3(B_grid2, T_grid2, G_grid_bt(:,:,5), [0 0], 'Color', [0 1 1], 'LineWidth', 2.5, 'DisplayName', 'Buckling Limit');
    % 
    % plot3(opt_b, opt_t, opt_W, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', 'DisplayName', 'Optimum');
    % plot3([opt_b, opt_b], [opt_t, opt_t], [0, opt_W], 'r:', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    % 
    % view(45, 20);
    % zlim([0, max(W_grid_bt(:))]);
    % xlabel('b (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('t (mm)', 'FontSize', 12, 'FontWeight', 'bold');
    % zlabel('W (kg)', 'FontSize', 12, 'FontWeight', 'bold');
    % legend('Location', 'eastoutside');
    % title(sprintf('b vs t Design Space\n(Fixed h=%.0f, \\theta=%.2f^\\circ)', opt_h, opt_th), 'FontSize', 14, 'FontWeight', 'bold');
    % grid on; hold off;
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function stop = outfun_fmincon(x, optimValues, state)
    % This function is called by fmincon at every iteration to log the path
    global hist_fmincon;
    stop = false;
    
    if isequal(state, 'iter') || isequal(state, 'init')
        % Format: [iteration, h, b, t, theta, objective_value]
        hist_fmincon = [hist_fmincon; optimValues.iteration, x(1), x(2), x(3), x(4), optimValues.fval];
    end
end