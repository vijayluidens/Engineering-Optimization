function [x_opt, f_opt, n_iter, history] = golden_section_search(f, a, b, tol)
% GOLDEN_SECTION_SEARCH  Self-implemented 1-D minimization
%
% Finds the minimum of f(x) on the interval [a, b] using the Golden
% Section Search method (based on Exercise 4 / Lecture 5 material).
%
% The interval is reduced by the golden ratio at each iteration until
% the bracket width falls below the tolerance.
%
% Inputs:
%   f    - function handle (scalar function of one variable)
%   a, b - initial bracket  [a < b]
%   tol  - convergence tolerance on bracket width
%
% Outputs:
%   x_opt   - location of minimum
%   f_opt   - function value at minimum
%   n_iter  - number of iterations taken
%   history - array of bracket widths at each iteration (for plotting)

    phi = (sqrt(5) - 1) / 2;   % golden ratio ~ 0.618
    
    % Initial interior points
    x1 = b - phi * (b - a);
    x2 = a + phi * (b - a);
    
    f1 = f(x1);
    f2 = f(x2);
    
    n_iter = 0;
    history = []; % Initialize array to store convergence data
    
    while (b - a) > tol
        n_iter = n_iter + 1;
        history(n_iter) = b - a; % Record the bracket width for this iteration
        
        if f1 < f2
            % Minimum is in [a, x2] — discard right section
            b  = x2;
            x2 = x1;   
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = f(x1);
        else
            % Minimum is in [x1, b] — discard left section
            a  = x1;
            x1 = x2;   
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = f(x2);
        end
    end
    
    x_opt = (a + b) / 2;
    f_opt = f(x_opt);
end