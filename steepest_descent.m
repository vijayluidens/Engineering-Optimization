function [x_opt, f_opt, history] = steepest_descent(obj, x0, lb, ub, options)
% STEEPEST_DESCENT  Self-implemented N-D steepest descent with line search
%
% Implements the steepest descent algorithm from Exercise 4:
%   1. Compute numerical gradient at current point  (finite differences)
%   2. Set search direction S = -gradient  (steepest descent)
%   3. Normalise S and compute alpha_max so search stays within bounds
%   4. Minimise along S using Golden Section Search (self-implemented)
%   5. Check convergence: |f_new - f_old| < TolF
%   6. Repeat
%
% Inputs:
%   obj     - objective function handle  f = obj(x)
%   x0      - starting point (column vector)
%   lb, ub  - lower/upper bounds on design variables
%   options - struct with fields:
%               TolF    : convergence tolerance on objective change
%               MaxIter : maximum number of iterations
%               h       : finite-difference step size (default 1e-7)
%
% Outputs:
%   x_opt   - optimal design variables found
%   f_opt   - objective value at optimum
%   history - [n_iter x (n+2)] matrix: [iter, x1, ..., xn, f] per row

    x   = x0(:);
    n   = length(x);
    h   = options.h;
    tol = options.TolF;

    history = zeros(options.MaxIter + 1, n + 2);
    f       = obj(x);
    history(1, :) = [0, x', f];

    for k = 1:options.MaxIter

        % ------------------------------------------------------------------
        % Step 1: Numerical gradient (forward finite differences)
        % ------------------------------------------------------------------
        grad = zeros(n, 1);
        for i = 1:n
            xe      = x;
            xe(i)   = xe(i) + h;
            grad(i) = (obj(xe) - f) / h;
        end

        % ------------------------------------------------------------------
        % Step 2: Search direction (steepest descent)
        % ------------------------------------------------------------------
        S      = -grad;
        S_norm = norm(S);
        if S_norm < 1e-14
            break;   % gradient is zero — already at a stationary point
        end
        S = S / S_norm;   % normalise to unit length

        % ------------------------------------------------------------------
        % Step 2.5: Search Direction Projection (Active Set Strategy)
        % ------------------------------------------------------------------
        % Prevent "Wall Paralysis" by zeroing out momentum pointing into bounds
        for i = 1:n
            % If at the upper bound and the search direction points higher
            if x(i) >= ub(i) - 1e-8 && S(i) > 0
                S(i) = 0.0;
                grad(i) = 0.0; % Clear gradient to prevent momentum buildup
            % If at the lower bound and the search direction points lower
            elseif x(i) <= lb(i) + 1e-8 && S(i) < 0
                S(i) = 0.0;
                grad(i) = 0.0; % Clear gradient to prevent momentum buildup
            end
        end

        % ------------------------------------------------------------------
        % Step 3: Compute alpha_max (largest step keeping x within bounds)
        % ------------------------------------------------------------------
        alpha_max = Inf;
        for i = 1:n
            if S(i) > 1e-15
                alpha_max = min(alpha_max, (ub(i) - x(i)) / S(i));
            elseif S(i) < -1e-15
                alpha_max = min(alpha_max, (lb(i) - x(i)) / S(i));
            end
        end
        
        if ~isfinite(alpha_max)
            alpha_max = 1.0;
        elseif alpha_max <= 0 
            alpha_max = 0; % Prevent ANY step if we are pushing into a wall
        end

        % ------------------------------------------------------------------
        % Step 4: Line search along S using Golden Section Search
        % ------------------------------------------------------------------
        line_obj    = @(alpha) obj(x + alpha * S);
        [alpha_opt, f_new] = golden_section_search(line_obj, 0, alpha_max, 1e-12);

        x_new = x + alpha_opt * S;
        x_new = max(lb, min(ub, x_new));

        % ------------------------------------------------------------------
        % Step 5: Store history
        % ------------------------------------------------------------------
        history(k + 1, :) = [k, x_new', f_new];

        % ------------------------------------------------------------------
        % Step 6: Convergence check
        % ------------------------------------------------------------------
        if abs(f_new - f) < tol
            history  = history(1:k+1, :);
            x        = x_new;
            f        = f_new;
            break;
        end

        x = x_new;
        f = f_new;
    end

    x_opt = x;
    f_opt = f;

end
