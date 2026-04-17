function phi = penalty_objective(x, p, mu)
% PENALTY_OBJECTIVE  Exterior penalty function for constrained optimization
%
%   phi(x, mu) = W(x) + mu * sum_i [ max(0, g_i(x)) ]^2
%              + mu * bound penalties
%
%   As mu -> infinity, the unconstrained minimum of phi converges to
%   the constrained minimum of W(x) subject to g_i(x) <= 0.
%
%   x = [h, b, t, theta]

    % Objective
    W = beam_weight(x, p);

    % Inequality constraints (g <= 0 is feasible)
    [g, ~] = beam_constraints(x, p);

    % Constraint violation penalty
    pen_con = 0;
    for i = 1:length(g)
        pen_con = pen_con + max(0, g(i))^2;
    end

    % Bound penalties — normalised by variable range
    pen_bnd = 0;
    for i = 1:length(x)
        range_i = p.ub(i) - p.lb(i);
        pen_bnd = pen_bnd + max(0, (p.lb(i) - x(i)) / range_i)^2;
        pen_bnd = pen_bnd + max(0, (x(i) - p.ub(i)) / range_i)^2;
    end

    phi = W + mu * (pen_con + pen_bnd);

end
