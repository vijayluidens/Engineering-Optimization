function W = beam_weight(x, p)
% BEAM_WEIGHT  Objective function: structural weight of the box beam spar (kg)
%
%   x = [h, b, t, theta]
%     h     = box height       (mm)
%     b     = box width        (mm)
%     t     = wall thickness   (mm)
%     theta = ply angle        (deg)  — does not affect weight
%
%   W = rho * perimeter * t * L = rho * 2*(h+b) * t * L

    h = x(1);
    b = x(2);
    t = x(3);
    % theta = x(4);   % not used in weight calculation

    W = p.rho * 2 * (h + b) * t * p.L;

end
