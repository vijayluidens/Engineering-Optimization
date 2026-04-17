function [x_opt, f_opt, history] = conjugate_gradient(obj, x0, lb, ub, options)
% CONJUGATE_GRADIENT  Self-implemented Conjugate Gradient descent
%
% Implements the Fletcher-Reeves Conjugate Gradient method as defined 
% in TU Delft Lecture 7, including mandatory N-step restarts to handle 
% non-quadratic functions and round-off errors.
%
% Inputs:
%   obj     - objective function handle  f = obj(x)
%   x0      - starting point (column vector)
%   lb, ub  - lower/upper bounds
%   options - struct with fields (TolF, MaxIter, h)

    x    = x0(:);
    n    = length(x); % n = 4 design variables
    h    = options.h;
    tol  = options.TolF;
    
    history = zeros(options.MaxIter + 1, n + 2);
    f       = obj(x);
    history(1, :) = [0, x', f];
    
    S_prev    = zeros(n, 1);   
    grad_prev = zeros(n, 1);   
    
    % NEW: Counter to track iterations since the last restart
    steps_since_restart = 0; 
    
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
        % Step 2: Conjugate direction & Mandatory Restarts
        % ------------------------------------------------------------------
        % Lecture 7, Slides 44 & 45: Restart every N steps to clear errors
        if k == 1 || steps_since_restart >= n
            beta = 0.0; % Pure steepest descent restart
            steps_since_restart = 0; % Reset counter
        else
            % BULLETPROOF FLETCHER-REEVES
            denom = grad_prev' * grad_prev;
            if denom < 1e-14
                beta = 0.0; % Safely restart if gradient was zeroed out
            else
                beta = (grad' * grad) / denom;
            end
        end
        
        % Search Direction
        S = -grad + beta * S_prev;
        
        % ------------------------------------------------------------------
        % Step 2.5: Search Direction Projection (Active Set Strategy)
        % ------------------------------------------------------------------
        % Prevent "Wall Paralysis" by zeroing out momentum pointing into bounds
        for i = 1:n
            if x(i) >= ub(i) - 1e-8 && S(i) > 0
                S(i) = 0.0;
                grad(i) = 0.0; 
            elseif x(i) <= lb(i) + 1e-8 && S(i) < 0
                S(i) = 0.0;
                grad(i) = 0.0; 
            end
        end
        
        % Safety: If conjugacy is entirely lost and points uphill, force immediate restart
        if dot(S, -grad) <= 0
            S = -grad;   
            steps_since_restart = 0;          
        end
        
        S_norm = norm(S);
        if S_norm < 1e-14
            break;   % stationary point reached
        end
        
        % ------------------------------------------------------------------
        % Step 3: Compute alpha_max using a NORMALIZED search direction
        % ------------------------------------------------------------------
        % Split-Personality: Normalize only for the line search to protect GSS precision
        S_search = S / S_norm;
        
        alpha_max = Inf;
        for i = 1:n
            if S_search(i) > 1e-15
                alpha_max = min(alpha_max, (ub(i) - x(i)) / S_search(i));
            elseif S_search(i) < -1e-15
                alpha_max = min(alpha_max, (lb(i) - x(i)) / S_search(i));
            end
        end
        
        if ~isfinite(alpha_max)
            alpha_max = 1.0;
        elseif alpha_max <= 0 
            alpha_max = 0; 
        end
        
        % ------------------------------------------------------------------
        % Step 4: Line search using Golden Section Search
        % ------------------------------------------------------------------
        line_obj          = @(alpha) obj(x + alpha * S_search);
        [alpha_opt, f_new] = golden_section_search(line_obj, 0, alpha_max, 1e-8);
        x_new             = x + alpha_opt * S_search;
        x_new             = max(lb, min(ub, x_new)); % Strict bounds clip
        
        % ------------------------------------------------------------------
        % Step 5: Step Rejection Fallback
        % ------------------------------------------------------------------
        if f_new >= f - 1e-12 % (Tiny buffer for floating point noise)
            S_sd = -grad;
            S_sd_norm = norm(S_sd);
            
            if S_sd_norm < 1e-14
                break;
            end
            
            S_sd_search = S_sd / S_sd_norm; % Normalize fallback for GSS
            
            % Recalculate boundary limits specifically for the new SD direction
            alpha_max_sd = Inf;
            for i = 1:n
                if S_sd_search(i) > 1e-15
                    alpha_max_sd = min(alpha_max_sd, (ub(i) - x(i)) / S_sd_search(i));
                elseif S_sd_search(i) < -1e-15
                    alpha_max_sd = min(alpha_max_sd, (lb(i) - x(i)) / S_sd_search(i));
                end
            end
            
            if ~isfinite(alpha_max_sd)
                alpha_max_sd = 1.0;
            elseif alpha_max_sd <= 0
                alpha_max_sd = 0; % Prevent wall breaches
            end
            
            line_obj_sd        = @(alpha) obj(x + alpha * S_sd_search);
            [alpha_opt, f_new] = golden_section_search(line_obj_sd, 0, alpha_max_sd, 1e-8);
            x_new              = x + alpha_opt * S_sd_search;
            x_new              = max(lb, min(ub, x_new)); % Strict bounds clip
            
            S                  = S_sd; % Keep the un-normalized vector for memory
            steps_since_restart = 0; % Reset counter due to forced restart
        end
        
        % ------------------------------------------------------------------
        % Step 6: Store history, prep for next iteration, check convergence
        % ------------------------------------------------------------------
        history(k + 1, :) = [k, x_new', f_new];
        if abs(f_new - f) < tol
            history = history(1:k+1, :);
            x       = x_new;
            f       = f_new;
            break;
        end
        
        S_prev    = S;      
        grad_prev = grad;   
        x         = x_new;
        f         = f_new;
        
        % Increment the restart counter!
        steps_since_restart = steps_since_restart + 1;
    end 
    
    x_opt = x;
    f_opt = f;
end