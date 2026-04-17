function [x_opt, f_opt, history] = slp_move_limits(obj, con, x0, lb, ub, ml, options)
% SLP_MOVE_LIMITS  Self-implemented Sequential Linear Programming with move limits
%
% Implements SLP from Exercise 7:
%   At each iteration around current point xk:
%     1. Linearise objective and constraints using finite-difference gradients
%     2. Solve the resulting LP with linprog (MATLAB built-in)
%     3. Update xk to the LP solution
%     4. Apply move limits: xk - ml <= x <= xk + ml  (intersected with lb/ub)
%     5. Check convergence: ||x_new - x_old|| < TolX
%
% Inputs:
%   obj     - objective function handle          f  = obj(x)
%   con     - constraint function handle         [c, ceq] = con(x)  (c <= 0)
%   x0      - starting point (column vector)
%   lb, ub  - design variable bounds
%   ml      - move limits vector (same size as x0)
%   options - struct with fields:
%               TolX    : convergence tolerance on ||dx||
%               MaxIter : maximum number of iterations
%               h       : finite-difference step size (default 1e-7)
%
% Outputs:
%   x_opt   - optimal design found
%   f_opt   - objective at optimum
%   history - [n_iter x (n+2)] matrix: [iter, x1, ..., xn, f] per row

    x   = x0(:);
    n   = length(x);
    h   = options.h;

    history    = zeros(options.MaxIter + 1, n + 2);
    f          = obj(x);
    history(1,:) = [0, x', f];

    for k = 1:options.MaxIter

        f_k = obj(x);

        % ------------------------------------------------------------------
        % Step 1: Compute gradients via forward finite differences
        % ------------------------------------------------------------------

        % Objective gradient
        df = zeros(n, 1);
        for i = 1:n
            xe    = x;  xe(i) = xe(i) + h;
            df(i) = (obj(xe) - f_k) / h;
        end

        % Constraint values and gradients
        [c_k, ~] = con(x);
        nc       = length(c_k);
        dc       = zeros(nc, n);    % each row = gradient of one constraint
        for j = 1:nc
            for i = 1:n
                xe         = x;  xe(i) = xe(i) + h;
                [c_e, ~]   = con(xe);
                dc(j, i)   = (c_e(j) - c_k(j)) / h;
            end
        end

        % ------------------------------------------------------------------
        % Step 2: Linearised problem
        %   min  df' * (x - xk)
        %   s.t. c_k + dc * (x - xk) <= 0
        %        lb_ml <= x <= ub_ml
        % ------------------------------------------------------------------

        f_lp = df;
        A_lp = dc;
        b_lp = -c_k(:) + dc * x;

        % Move-limited bounds
        lb_ml = max(lb(:), x - ml(:));
        ub_ml = min(ub(:), x + ml(:));

        opts_lp = optimoptions('linprog', 'Display', 'off');

        try
            [x_new, ~, exitflag] = linprog(f_lp, A_lp, b_lp, [], [], lb_ml, ub_ml, opts_lp);
        catch
            exitflag = -1;
        end

        if isempty(x_new) || exitflag < 0
            x_new = x;   % stay at current point if LP infeasible
        end

        % ------------------------------------------------------------------
        % Step 3: Evaluate new objective and store
        % ------------------------------------------------------------------
        f_new = obj(x_new);
        history(k+1, :) = [k, x_new', f_new];

        % ------------------------------------------------------------------
        % Step 4: Convergence check
        % ------------------------------------------------------------------
        if norm(x_new - x) < options.TolX
            history = history(1:k+1, :);
            x = x_new;
            break;
        end

        x = x_new;
    end

    x_opt = x;
    f_opt = obj(x);

end
