function [Ex, Ey, Gxy, nu_xy] = compute_effective_properties(theta_deg, p)
% COMPUTE_EFFECTIVE_PROPERTIES  Classical Lamination Theory for [0/theta/-theta/0]s
%
% Computes effective in-plane engineering constants of the symmetric
% balanced laminate using the A-matrix (extensional stiffness).
%
% Inputs:
%   theta_deg — off-axis ply angle in degrees
%   p         — parameter struct (must contain E1, E2, G12, nu12, t_ply)
%
% Outputs:
%   Ex, Ey  — effective axial moduli        (MPa)
%   Gxy     — effective in-plane shear modulus (MPa)
%   nu_xy   — effective major Poisson's ratio

    E1  = p.E1;   E2  = p.E2;
    G12 = p.G12;  nu12 = p.nu12;
    nu21 = nu12 * E2 / E1;

    % On-axis reduced stiffness matrix
    denom = 1 - nu12 * nu21;
    Q11 = E1  / denom;
    Q22 = E2  / denom;
    Q12 = nu12 * E2 / denom;
    Q66 = G12;

    % Ply angles: [0/theta/-theta/0]s = [0, th, -th, 0, 0, -th, th, 0]
    angles = [0, theta_deg, -theta_deg, 0, 0, -theta_deg, theta_deg, 0];

    % Build A matrix using explicit Q-bar transformation formulas
    A11 = 0;  A22 = 0;  A12 = 0;  A66 = 0;

    for k = 1:length(angles)
        th = angles(k) * pi / 180;
        c  = cos(th);  s = sin(th);
        c2 = c^2;  s2 = s^2;  c4 = c^4;  s4 = s^4;  c2s2 = c2 * s2;

        % Transformed reduced stiffness (Tsai convention)
        Qb11 = Q11*c4 + 2*(Q12 + 2*Q66)*c2s2 + Q22*s4;
        Qb22 = Q11*s4 + 2*(Q12 + 2*Q66)*c2s2 + Q22*c4;
        Qb12 = (Q11 + Q22 - 4*Q66)*c2s2 + Q12*(c4 + s4);
        Qb66 = (Q11 + Q22 - 2*Q12 - 2*Q66)*c2s2 + Q66*(c4 + s4);

        A11 = A11 + Qb11 * p.t_ply;
        A22 = A22 + Qb22 * p.t_ply;
        A12 = A12 + Qb12 * p.t_ply;
        A66 = A66 + Qb66 * p.t_ply;
    end

    % Total laminate thickness
    t_lam = length(angles) * p.t_ply;

    % Compliance (invert 3x3 A-matrix; A16 = A26 = 0 for balanced layup)
    A_mat = [A11, A12, 0;
             A12, A22, 0;
             0,   0,   A66];
    a = inv(A_mat);

    % Effective engineering constants
    Ex   = 1 / (a(1,1) * t_lam);
    Ey   = 1 / (a(2,2) * t_lam);
    Gxy  = 1 / (a(3,3) * t_lam);
    nu_xy = -a(1,2) / a(1,1);

end
