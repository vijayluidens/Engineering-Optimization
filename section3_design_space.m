function section3_design_space(p)
% SECTION3_DESIGN_SPACE  Parameter interaction analysis for 4-variable problem
%
% Since the design space is 4-D, we show 2-D slices

    N = 200;

    % =====================================================================
    % Slice 1: h vs b  (fixed t = 0.5 mm, theta = 45 deg)
    % =====================================================================
    t_fix  = 1.5;
    th_fix = 60;

    h_vec = linspace(p.lb(1), p.ub(1), N);
    b_vec = linspace(p.lb(2), p.ub(2), N);
    [H, B] = meshgrid(h_vec, b_vec);

    W_hb  = p.rho * 2 * (H + B) * t_fix * p.L;

    % Evaluate all 5 constraints on grid
    G_hb = zeros(N, N, 5);
    for ii = 1:N
        for jj = 1:N
            [gg, ~] = beam_constraints([H(ii,jj), B(ii,jj), t_fix, th_fix], p);
            G_hb(ii,jj,:) = gg;
        end
    end

    % -----------------------------------------------------------------------
    % Figure 1: Objective contours W(h, b) — 2D contour with ShowText
    % -----------------------------------------------------------------------
    figure('Name', 'Fig1: Weight Contours h-b', 'Position', [50 50 700 500]);
    contour(H, B, W_hb, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;
    % Constraint boundaries (assignment color scheme: r, g, k, b, c)
    contour(H, B, G_hb(:,:,1), [0 0], 'r', 'LineWidth', 1.5);
    contour(H, B, G_hb(:,:,3), [0 0], 'k', 'LineWidth', 1.5);
    contour(H, B, G_hb(:,:,4), [0 0], 'b', 'LineWidth', 1.5);
    contour(H, B, G_hb(:,:,5), [0 0], 'c', 'LineWidth', 1.5);
    xlabel('Height h (mm)');
    ylabel('Width b (mm)');
    title(sprintf('W(h, b) at t = %.1f mm, \\theta = %.0f deg', t_fix, th_fix));
    legend({'W (kg)', 'g_1: Bending stress', 'g_3: Tip deflection', ...
            'g_4: Twist angle', 'g_5: Web buckling'}, 'Location', 'eastoutside');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 2: Design space (h vs b) — contour + constraint boundaries
    % -----------------------------------------------------------------------
    
    figure('Name', 'Fig2: Design Space h-b', 'Position', [100 50 750 550]);

    % Objective contours with text labels
    contour(H, B, W_hb, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;

    % Constraint boundaries at g=0 (solid, LineWidth 1.5)
    contour(H, B, G_hb(:,:,1), [0 0], 'r', 'LineWidth', 1.5, ...
            'DisplayName', 'g_1: stress');
    contour(H, B, G_hb(:,:,3), [0 0], 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'g_3: deflection');
    contour(H, B, G_hb(:,:,4), [0 0], 'b', 'LineWidth', 1.5, ...
            'DisplayName', 'g_4: twist');
    contour(H, B, G_hb(:,:,5), [0 0], 'c', 'LineWidth', 1.5, ...
            'DisplayName', 'g_5: buckling');

    % Dashed infeasible-side indicators (slightly inside violation region)
    contour(H, B, G_hb(:,:,1), [0.05 0.05], '--r', 'LineWidth', 1);
    contour(H, B, G_hb(:,:,3), [0.05 0.05], '--k', 'LineWidth', 1);
    contour(H, B, G_hb(:,:,4), [0.05 0.05], '--b', 'LineWidth', 1);
    contour(H, B, G_hb(:,:,5), [0.05 0.05], '--c', 'LineWidth', 1);

    xlabel('Height h (mm)');
    ylabel('Width b (mm)');
    title(sprintf('Design space (h, b) at t = %.1f mm, \\theta = %.0f deg', t_fix, th_fix));
    legend({'W (kg)', 'g_1: stress', 'g_3: deflection', ...
            'g_4: twist', 'g_5: buckling'}, 'Location', 'northwest');
    grid on; hold off;

    % =====================================================================
    % Slice 2: t vs theta  (fixed h = 50 mm, b = 30 mm)
    % =====================================================================
    h_fix = 50;
    b_fix = 30;

    t_vec  = linspace(p.lb(3), p.ub(3), N);
    th_vec = linspace(p.lb(4)+0.5, p.ub(4)-0.5, N);   % avoid exact 0 and 90
    [Tg, TH] = meshgrid(t_vec, th_vec);

    W_tth = p.rho * 2 * (h_fix + b_fix) * Tg * p.L;

    G_tth = zeros(N, N, 5);
    for ii = 1:N
        for jj = 1:N
            [gg, ~] = beam_constraints([h_fix, b_fix, Tg(ii,jj), TH(ii,jj)], p);
            G_tth(ii,jj,:) = gg;
        end
    end

    % -----------------------------------------------------------------------
    % Figure 3: Objective contours W(t, theta)
    % -----------------------------------------------------------------------
    figure('Name', 'Fig3: Weight Contours t-theta', 'Position', [150 50 700 500]);
    contour(Tg, TH, W_tth, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;
    contour(Tg, TH, G_tth(:,:,1), [0 0], 'r', 'LineWidth', 1.5);
    contour(Tg, TH, G_tth(:,:,2), [0 0], 'g', 'LineWidth', 1.5);
    contour(Tg, TH, G_tth(:,:,3), [0 0], 'k', 'LineWidth', 1.5);
    contour(Tg, TH, G_tth(:,:,4), [0 0], 'b', 'LineWidth', 1.5);
    contour(Tg, TH, G_tth(:,:,5), [0 0], 'c', 'LineWidth', 1.5);
    xlabel('Thickness t (mm)');
    ylabel('Ply angle \theta (deg)');
    title(sprintf('W(t, \\theta) at h = %.0f mm, b = %.0f mm', h_fix, b_fix));
    legend({'W (kg)', 'g_1: stress', 'g_2: shear', 'g_3: deflection', ...
            'g_4: twist', 'g_5: buckling'}, 'Location', 'northeast');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 4: Design space (t vs theta)
    % -----------------------------------------------------------------------
    figure('Name', 'Fig4: Design Space t-theta', 'Position', [200 50 750 550]);

    contour(Tg, TH, W_tth, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;

    contour(Tg, TH, G_tth(:,:,1), [0 0], 'r', 'LineWidth', 1.5, ...
            'DisplayName', 'g_1: stress');
    contour(Tg, TH, G_tth(:,:,2), [0 0], 'g', 'LineWidth', 1.5, ...
            'DisplayName', 'g_2: shear');
    contour(Tg, TH, G_tth(:,:,3), [0 0], 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'g_3: deflection');
    contour(Tg, TH, G_tth(:,:,4), [0 0], 'b', 'LineWidth', 1.5, ...
            'DisplayName', 'g_4: twist');
    contour(Tg, TH, G_tth(:,:,5), [0 0], 'c', 'LineWidth', 1.5, ...
            'DisplayName', 'g_5: buckling');

    % Dashed infeasible-side indicators
    contour(Tg, TH, G_tth(:,:,1), [0.05 0.05], '--r', 'LineWidth', 1);
    contour(Tg, TH, G_tth(:,:,2), [0.05 0.05], '--g', 'LineWidth', 1);
    contour(Tg, TH, G_tth(:,:,3), [0.05 0.05], '--k', 'LineWidth', 1);
    contour(Tg, TH, G_tth(:,:,4), [0.05 0.05], '--b', 'LineWidth', 1);
    contour(Tg, TH, G_tth(:,:,5), [0.05 0.05], '--c', 'LineWidth', 1);

    xlabel('Thickness t (mm)');
    ylabel('Ply angle \theta (deg)');
    title(sprintf('Design space (t, \\theta) at h = %.0f mm, b = %.0f mm', h_fix, b_fix));
    legend({'W (kg)', 'g_1: stress', 'g_2: shear', 'g_3: deflection', ...
            'g_4: twist', 'g_5: buckling'}, 'Location', 'northwest');
    grid on; hold off;

    % =====================================================================
    % Slice 3: h vs t  (fixed b = 30 mm, theta = 60 deg)
    % =====================================================================
    figure('Name', 'Fig5: Design Space h-t', 'Position', [250 50 750 550]);

    [Hg2, Tg2] = meshgrid(h_vec, t_vec);
    W_ht = p.rho * 2 * (Hg2 + b_fix) .* Tg2 * p.L;

    G_ht = zeros(N, N, 5);
    for ii = 1:N
        for jj = 1:N
            [gg, ~] = beam_constraints([Hg2(ii,jj), b_fix, Tg2(ii,jj), th_fix], p);
            G_ht(ii,jj,:) = gg;
        end
    end

    contour(Hg2, Tg2, W_ht, 15, 'ShowText', 'on', 'LabelSpacing', 500);
    hold on;

    contour(Hg2, Tg2, G_ht(:,:,1), [0 0], 'r', 'LineWidth', 1.5, ...
            'DisplayName', 'g_1: stress');
    contour(Hg2, Tg2, G_ht(:,:,3), [0 0], 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'g_3: deflection');
    contour(Hg2, Tg2, G_ht(:,:,4), [0 0], 'b', 'LineWidth', 1.5, ...
            'DisplayName', 'g_4: twist');
    contour(Hg2, Tg2, G_ht(:,:,5), [0 0], 'c', 'LineWidth', 1.5, ...
            'DisplayName', 'g_5: buckling');

    % Dashed infeasible-side indicators
    contour(Hg2, Tg2, G_ht(:,:,1), [0.05 0.05], '--r', 'LineWidth', 1);
    contour(Hg2, Tg2, G_ht(:,:,3), [0.05 0.05], '--k', 'LineWidth', 1);
    contour(Hg2, Tg2, G_ht(:,:,4), [0.05 0.05], '--b', 'LineWidth', 1);
    contour(Hg2, Tg2, G_ht(:,:,5), [0.05 0.05], '--c', 'LineWidth', 1);

    xlabel('Height h (mm)');
    ylabel('Thickness t (mm)');
    title(sprintf('Design space (h, t) at b = %.0f mm, \\theta = %.0f deg', b_fix, th_fix));
    legend({'W (kg)', 'g_1: stress', 'g_3: deflection', ...
            'g_4: twist', 'g_5: buckling'}, 'Location', 'northwest');
    grid on; hold off;

    % -----------------------------------------------------------------------
    % Figure 6: Sensitivity maps
    % -----------------------------------------------------------------------
    figure('Name', 'Fig6: Logarithmic Sensitivity', 'Position', [100 50 1500 420]);

    h_vec = linspace(p.lb(1), p.ub(1), 200);
    b_vec = linspace(p.lb(2), p.ub(2), 200);
    [H_g, B_g] = meshgrid(h_vec, b_vec);
    
    SL_h  =  H_g ./ (H_g + B_g);
    SL_b  =  B_g ./ (H_g + B_g);
    SL_diff = (H_g - B_g) ./ (H_g + B_g);   % dominance map
    
    % ── Panel 1: S_L(h) ───────────────────────────────────────────────────────
    subplot(1,5,1);
    contourf(H_g, B_g, SL_h, 20);
    colorbar; clim([0 1]);
    xlabel('h (mm)', 'FontWeight', 'bold'); ylabel('b (mm)', 'FontWeight', 'bold');
    title({'S_L(h) = h/(h+b)'; 'h-b design space'});
    grid on;
    
    % ── Panel 2: S_L(b) ───────────────────────────────────────────────────────
    subplot(1,5,2);
    contourf(H_g, B_g, SL_b, 20);
    colorbar; clim([0 1]);
    xlabel('h (mm)', 'FontWeight', 'bold'); ylabel('b (mm)', 'FontWeight', 'bold');
    title({'S_L(b) = b/(h+b)'; 'h-b design space'});
    grid on;
    
    % ── Panel 3: Dominance map S_L(h) - S_L(b) ───────────────────────────────
    subplot(1,5,3);
    contourf(H_g, B_g, SL_diff, 20);
    colorbar; clim([-1 1]);
    hold on;
    contour(H_g, B_g, SL_diff, [0 0], 'w-', 'LineWidth', 2.5);  % h=b line
    n = 256;
    cmap = [linspace(0.8,1,n/2)', linspace(0.1,1,n/2)', linspace(0.1,1,n/2)'; ...
        linspace(1,0.1,n/2)', linspace(1,0.1,n/2)', linspace(1,0.8,n/2)'];
    colormap(gca, cmap);  % red=h dominates, blue=b dominates
    xlabel('h (mm)', 'FontWeight', 'bold'); ylabel('b (mm)', 'FontWeight', 'bold');
    title({'S_L(h) - S_L(b) = (h-b)/(h+b)'; ''}); 
    % text(mean(h_vec)*1.1, mean(b_vec)*1.35, 'h dominant', ...
    %      'Color','w','FontWeight','bold','FontSize',8,'HorizontalAlignment','center');
    % text(mean(h_vec)*0.65, mean(b_vec)*0.75, 'b dominant', ...
    %      'Color','w','FontWeight','bold','FontSize',8,'HorizontalAlignment','center');
    grid on; hold off;
    
    % ── Panel 4: Bar chart at nominal ─────────────────────────────────────────
    subplot(1,5,4);
    h0 = mean([p.lb(1), p.ub(1)]); 
    b0 = mean([p.lb(2), p.ub(2)]); 
    SL_vals = [h0/(h0+b0), b0/(h0+b0), 1, 0];
    labels  = {'h', 'b', 't', '\theta'};
    colors  = [0 0.45 0.74; 0.85 0.33 0.1; 0.47 0.67 0.19; 0.7 0.7 0.7];
    b_hdl = bar(SL_vals, 'FaceColor', 'flat');
    b_hdl.CData = colors;
    set(gca, 'XTickLabel', labels, 'FontSize', 10, 'FontWeight', 'bold');
    yline(1, 'r--', 'LineWidth', 1.5);
    ylabel('S_L  [ - ]', 'FontWeight', 'bold');
    title({'Nominal Bar Chart'; sprintf('h=%.0f, b=%.0f mm', h0, b0)}, 'FontSize', 10);
    ylim([0 1.3]);
    for k = 1:4
        text(k, SL_vals(k)+0.04, sprintf('%.3f', SL_vals(k)), ...
             'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
    end
    % text(2.5, 1.12, 'S_L(t)=1: unit influence', 'FontSize',8,'Color','r',...
    %      'HorizontalAlignment','center');
    grid on;
    
    % ── Panel 5: Trade-off along h sweep (b fixed at nominal) ─────────────────
    subplot(1,5,5);
    SL_h_line = h_vec ./ (h_vec + b0);
    SL_b_line = b0    ./ (h_vec + b0);
    plot(h_vec, SL_h_line, 'b-',  'LineWidth', 2, 'DisplayName', 'S_L(h)');
    hold on;
    plot(h_vec, SL_b_line, 'r-', 'LineWidth', 2, 'DisplayName', 'S_L(b)');
    yline(1, 'k:',  'LineWidth', 1.5, 'DisplayName', 'S_L(t) = 1');
    yline(0, 'm:',  'LineWidth', 1.5, 'DisplayName', 'S_L(\theta) = 0');
    xline(h0, 'k-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xline(b0, 'k--','LineWidth', 1.0, 'HandleVisibility', 'off');  % h=b0 crossover
    xlabel('h (mm)', 'FontWeight', 'bold'); ylabel('S_L  [ - ]', 'FontWeight', 'bold');
    title({'S_L Trade-off vs. h'});
    legend('Location', 'east', 'FontSize', 9);
    ylim([-0.1 1.3]); 
    xlim([p.lb(1), p.ub(1)]);
    grid on; hold off;
    
    % sgtitle({'Logarithmic Sensitivity Analysis:   S_L = (x_i / W) \cdot \partialW/\partialx_i'; ...
    %          'Dimensionless — directly comparable across all variables'}, ...
    %         'FontSize', 12, 'FontWeight', 'bold');
    fprintf('  Figures 1-6 generated.\n');
end
