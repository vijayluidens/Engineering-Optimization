function section_monotonicity(p)
% SECTION_MONOTONICITY  Visualise that W is strictly increasing in h, b, t
%
% Layout (2x3):
%   Row 1 — 2-D contour slices (pairs of variables, third fixed at nominal)
%   Row 2 — 1-D slices (one variable swept, other two fixed at nominal)

    % ── Nominal (fixed) values ─────────────────────────────────────────────
    h0  = 50;    % mm  – mid-range nominal
    b0  = 30;    % mm
    t0  = 2;     % mm
    th0 = 45;    % deg (theta only enters constraints, not W)

    % ── Sweep vectors ──────────────────────────────────────────────────────
    h_vec = linspace(p.lb(1), p.ub(1), 200);
    b_vec = linspace(p.lb(2), p.ub(2), 200);
    t_vec = linspace(p.lb(3), p.ub(3), 200);

    % ── Objective formula  W = rho * 2*(h+b) * t * L ──────────────────────
    W_h  = p.rho * 2*(h_vec + b0) .* t0 .* p.L;   % vary h
    W_b  = p.rho * 2*(h0  + b_vec).* t0 .* p.L;   % vary b
    W_t  = p.rho * 2*(h0  + b0)   .* t_vec .* p.L; % vary t

    % 2-D grids
    [H,  B]  = meshgrid(h_vec, b_vec);
    [H2, T2] = meshgrid(h_vec, t_vec);
    [B2, T3] = meshgrid(b_vec, t_vec);

    W_hb = p.rho * 2*(H  + B)  .* t0 .* p.L;
    W_ht = p.rho * 2*(H2 + b0) .* T2 .* p.L;
    W_bt = p.rho * 2*(h0 + B2) .* T3 .* p.L;

    % ======================================================================
    figure('Name', 'Monotonicity Analysis', 'Position', [50 50 1200 700]);

    % ── Row 1: 2-D contours ───────────────────────────────────────────────
    subplot(2,3,1);
    contour(H, B, W_hb, 15, 'ShowText', 'on', 'LabelSpacing', 400);
    xlabel('h (mm)'); ylabel('b (mm)');
    title(sprintf('W(h,b),  t = %.0f mm', t0));
    colorbar; grid on;

    subplot(2,3,2);
    contour(H2, T2, W_ht, 15, 'ShowText', 'on', 'LabelSpacing', 400);
    xlabel('h (mm)'); ylabel('t (mm)');
    title(sprintf('W(h,t),  b = %.0f mm', b0));
    colorbar; grid on;

    subplot(2,3,3);
    contour(B2, T3, W_bt, 15, 'ShowText', 'on', 'LabelSpacing', 400);
    xlabel('b (mm)'); ylabel('t (mm)');
    title(sprintf('W(b,t),  h = %.0f mm', h0));
    colorbar; grid on;

    % ── Row 2: 1-D slices — the core monotonicity evidence ────────────────
    ax4 = subplot(2,3,4);
    plot(h_vec, W_h, 'b-', 'LineWidth', 2);
    xlabel('h (mm)'); ylabel('W (kg)');
    title(sprintf('W vs h  |  b=%.0f, t=%.0f', b0, t0));
    annotation_arrow(ax4, h_vec, W_h, '\partialW/\partialh > 0');
    grid on;

    ax5 = subplot(2,3,5);
    plot(b_vec, W_b, 'r-', 'LineWidth', 2);
    xlabel('b (mm)'); ylabel('W (kg)');
    title(sprintf('W vs b  |  h=%.0f, t=%.0f', h0, t0));
    annotation_arrow(ax5, b_vec, W_b, '\partialW/\partialb > 0');
    grid on;

    ax6 = subplot(2,3,6);
    plot(t_vec, W_t, 'k-', 'LineWidth', 2);
    xlabel('t (mm)'); ylabel('W (kg)');
    title(sprintf('W vs t  |  h=%.0f, b=%.0f', h0, b0));
    annotation_arrow(ax6, t_vec, W_t, '\partialW/\partialt > 0');
    grid on;

    sgtitle('Monotonicity Analysis: W strictly increases in h, b, and t', ...
            'FontSize', 13, 'FontWeight', 'bold');

    fprintf('  Monotonicity figure generated.\n');
end


% ── Helper: add a small upward arrow + label in the top-left of an axes ──
function annotation_arrow(ax, ~, ~, label_str)
    % Place a text label inside the axes
    xl = xlim(ax);  yl = ylim(ax);
    text(ax, xl(1) + 0.05*(xl(2)-xl(1)), ...
             yl(2) - 0.08*(yl(2)-yl(1)), ...
             label_str, ...
             'FontSize', 11, 'Color', [0.2 0.2 0.8], ...
             'FontWeight', 'bold');
end