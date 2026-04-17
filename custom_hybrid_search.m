function [x_opt, f_opt, iter, history] = custom_hybrid_search(f, a, b, tol)
    % CUSTOM_HYBRID_SEARCH: A self-made 1D optimization algorithm.
    % Stage 1: Golden Section Search for reliable bracketing.
    % Stage 2: Successive Parabolic Interpolation for rapid local convergence.

    phi = (sqrt(5) - 1) / 2;
    iter = 0;
    max_iter = 100;
    history = []; % Array to store the bracket width/error for plotting

    % ==========================================
    % STAGE 1: Golden Section Search (GSS)
    % ==========================================
    initial_width = b - a;
    x1 = b - phi * (b - a);
    x2 = a + phi * (b - a);
    f1 = f(x1);
    f2 = f(x2);

    % Run GSS until the interval shrinks by a factor of 100
    while (b - a) > 0.01 * initial_width && iter < max_iter
        iter = iter + 1;
        history(iter) = b - a; % Record interval width

        if f1 < f2
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = f(x1);
        else
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = f(x2);
        end
    end

    % ==========================================
    % STAGE 2: Successive Parabolic Interpolation
    % ==========================================
    % Identify the 3 best points from our final GSS bracket
    points = [a, x1, x2, b];
    f_vals = [f(a), f1, f2, f(b)];
    [~, sort_idx] = sort(f_vals);

    % Assign x1, x2, x3 based on best function values (x1 is the lowest)
    x_1 = points(sort_idx(1)); f_1 = f_vals(sort_idx(1));
    x_2 = points(sort_idx(2)); f_2 = f_vals(sort_idx(2));
    x_3 = points(sort_idx(3)); f_3 = f_vals(sort_idx(3));

    current_width = max([x_1, x_2, x_3]) - min([x_1, x_2, x_3]);

    % Run SPI until the maximum distance between our 3 points is below tolerance
    while current_width > tol && iter < max_iter
        iter = iter + 1;
        history(iter) = current_width;

        % Parabolic Interpolation Formula
        num = (x_2 - x_1)^2 * (f_2 - f_3) - (x_2 - x_3)^2 * (f_2 - f_1);
        den = (x_2 - x_1) * (f_2 - f_3) - (x_2 - x_3) * (f_2 - f_1);

        % Guard against division by zero (collinear points)
        if abs(den) < 1e-14
            break; 
        end

        % Calculate new point
        x_p = x_2 - 0.5 * (num / den);
        
        % Prevent algorithm stalling if step size is virtually zero
        if abs(x_p - x_1) < 1e-10
            break;
        end
        
        f_p = f(x_p);

        % Combine old points and new point, then keep the best 3
        all_x = [x_1, x_2, x_3, x_p];
        all_f = [f_1, f_2, f_3, f_p];
        [~, s_idx] = sort(all_f);

        x_1 = all_x(s_idx(1)); f_1 = all_f(s_idx(1));
        x_2 = all_x(s_idx(2)); f_2 = all_f(s_idx(2));
        x_3 = all_x(s_idx(3)); f_3 = all_f(s_idx(3));

        current_width = max([x_1, x_2, x_3]) - min([x_1, x_2, x_3]);
    end

    % Return the best point found
    x_opt = x_1;
    f_opt = f_1;
end