% RUN_ALL.M  3D Composite Box Beam — Structural Optimization
%
% Mirrors the report structure:
%   Initial Investigation - Boundedness, monotonicity, convexity, noise
%   Section 3 - Design space / parameter interaction analysis (2D slices)
%   Section 4 - Simplified optimization (2D and 1D reductions)
%   Section 5 - Full 4-D optimization (fmincon, penalty, SD, SLP)
%
% Model: thin-walled rectangular box beam with CFRP composite walls
%   Design variables: x = [h, b, t, theta]
%   Constraints: bending stress, torsional shear, deflection, twist, buckling
%
% Run this script from the 'vijay model' folder, or add it to the path.
% All figures are saved automatically to the 'figures/' subfolder.
%
% =========================================================================

clear; clc; close all;

% Add this folder to path so helper functions are visible
addpath(fileparts(mfilename('fullpath')));

% Create figures output folder
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% -------------------------------------------------------------------------
% Load parameters
% -------------------------------------------------------------------------
p = get_parameters();

fprintf('=============================================\n');
fprintf(' 3D Composite Box Beam Optimization\n');
fprintf(' (F1 Rear Wing Spar)\n');
fprintf('=============================================\n\n');
fprintf('Design variables: x = [h, b, t, theta]\n');
fprintf('  h     = spar box height      [%.0f, %.0f] mm\n', p.lb(1), p.ub(1));
fprintf('  b     = spar box width       [%.0f, %.0f] mm\n', p.lb(2), p.ub(2));
fprintf('  t     = wall thickness       [%.1f, %.1f] mm\n', p.lb(3), p.ub(3));
fprintf('  theta = off-axis ply angle   [%.0f, %.0f] deg\n', p.lb(4), p.ub(4));
fprintf('\nFixed parameters:\n');
fprintf('  Span length        L = %.0f mm\n',       p.L);
fprintf('  Tip vertical load  P = %.0f N\n',        p.P);
fprintf('  Tip torque         T = %.0f N*mm (%.0f Nm)\n', p.T_torque, p.T_torque/1e3);
fprintf('  CFRP density       rho = %.0f kg/m^3\n', p.rho/1e-9);
fprintf('  Laminate           [0/theta/-theta/0]s,  %d plies x %.3f mm\n', p.n_plies, p.t_ply);
fprintf('  E1 = %.0f GPa,  E2 = %.1f GPa,  G12 = %.2f GPa,  nu12 = %.2f\n', ...
        p.E1/1e3, p.E2/1e3, p.G12/1e3, p.nu12);
fprintf('\nStructural allowables:\n');
fprintf('  sigma_allow = %.0f MPa\n',  p.sigma_allow);
fprintf('  tau_allow   = %.0f MPa\n',  p.tau_allow);
fprintf('  delta_max   = %.0f mm\n',   p.delta_max);
fprintf('  phi_max     = %.1f deg\n',  p.phi_max*180/pi);
fprintf('\n');

% -------------------------------------------------------------------------
% Initial Investigation (boundedness, monotonicity, convexity, noise)
% -------------------------------------------------------------------------
fprintf('-----------------------------------------\n');
fprintf(' Initial Problem Investigation\n');
fprintf('-----------------------------------------\n');
section_investigation(p);

% -------------------------------------------------------------------------
% Section 3: Design space analysis
% -------------------------------------------------------------------------
fprintf('\n-----------------------------------------\n');
fprintf(' Section 3: Design Space Analysis\n');
fprintf('-----------------------------------------\n');
section3_design_space(p);

% -------------------------------------------------------------------------
% Section 4: Simplified optimization
% -------------------------------------------------------------------------
fprintf('\n-----------------------------------------\n');
fprintf(' Section 4: Simplified Optimization\n');
fprintf('-----------------------------------------\n');
section4_simplified(p);

% -------------------------------------------------------------------------
% Section 5: Full 4-D constrained optimization
% -------------------------------------------------------------------------
fprintf('\n-----------------------------------------\n');
fprintf(' Section 5: Full 4-D Optimization\n');
fprintf('-----------------------------------------\n');
section5_optimization(p);

% -------------------------------------------------------------------------
% Save all figures (PNG for screen + PDF for report, matching assignment style)
% -------------------------------------------------------------------------
fprintf('\nSaving figures to:  %s\n', fig_dir);
fig_handles = findall(groot, 'Type', 'figure');
for i = 1:length(fig_handles)
    fname = fig_handles(i).Name;
    fname = regexprep(fname, '[^a-zA-Z0-9_]', '_');
    if isempty(fname), fname = sprintf('figure_%d', i); end
    saveas(fig_handles(i), fullfile(fig_dir, [fname '.png']));
    saveas(fig_handles(i), fullfile(fig_dir, [fname '.pdf']));
end
fprintf('Saved %d figures (PNG + PDF).\n', length(fig_handles));

fprintf('\n=============================================\n');
fprintf(' Done.\n');
fprintf('=============================================\n');
