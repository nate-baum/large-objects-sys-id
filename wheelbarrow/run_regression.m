clear;
clc;
close all;

addpath('../utils/plot');
addpath('../utils/log');
addpath('../utils/rigid_body');
load('results/wheelbarrow_data.mat', 'data', 'phi_sw', 'dt');
log_msg('run_regression', 'loaded wheelbarrow_data.mat (%d samples)', data.N);

d1 = 1879:4883;
d2 = 6000:7755;
d3 = 9265:11646;
d4 = 12422:14345;
d5 = 15713:17817;
d6 = 18586:20859;

indices_train   = d1;
indices_test    = d6;


m_known = 45.0;
lambda  = 100;
mu_r    = 0.02;

sigma_f = [950 * 0.015; 950 * 0.015; 1900 * 0.010];
sigma_n = [40 * 0.0125; 40 * 0.0125; 40 * 0.0125];

data.C = [zeros(3, 3), eye(3)];
data.G = diag(sigma_n(1) ./ [sigma_f; sigma_n]);

data.b = [-13.7; 0; 0; 0; 0; 9.5];

training_data = get_data_selection(data, indices_train);
testing_data = get_data_selection(data, indices_test);
log_msg('run_regression', 'training set: %d samples, testing set: %d samples', ...
    training_data.N, testing_data.N);

theta_CAD = phi_to_logchol(phi_sw);
x0 = [theta_CAD(2:10); log(0.5); log(0.5); log(1.0)];

opts = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 10000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'Display', 'iter');

cost = @(x) evaluate_cost_function(x, training_data, m_known, x0, lambda);
log_msg('run_regression', 'starting optimization...');
[x_opt, resnorm] = lsqnonlin(cost, x0, [], [], opts);
log_msg('run_regression', 'optimization finished, resnorm=%.4g', resnorm);

[phi_final, bv, bc, bpsi] = extract_wheelbarrow_results(x_opt, m_known);

J_sol = compute_jacobian(x_opt, training_data, m_known, 0.0);

[~, S, ~] = svd(J_sol, 'econ');
s = diag(S);
kappa = cond(S);

if kappa > 1e6
    log_msg('run_regression', 'warning: cond(J_data)=%.2e is poorly conditioned', kappa);
end

kappa_eff = sqrt(s(1) ^ 2 + lambda ^ 2) / sqrt(s(end) ^ 2 + lambda ^ 2);

% -------------------------------------------------------------------------
% Compute test-set predictions (needed for both the report table and the plot)
% -------------------------------------------------------------------------
N_all = data.N;
t_col_all = (0:N_all - 1) * dt;

N_test = testing_data.N;
w_hand = zeros(6, N_test);
w_bal = zeros(6, N_test);

for k = 1:N_test

    f = evaluate_f(testing_data, k, phi_final, bv, bc, bpsi);

    w_hand(:, k) = testing_data.y(:, k);
    w_bal(:, k) = f + testing_data.b;

end

resid = w_hand - w_bal;

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
log_msg('run_regression', 'Training set: d%d   Test set: d%d', ...
    find(cellfun(@(d) isequal(d, indices_train), {d1, d2, d3, d4, d5, d6})), ...
    find(cellfun(@(d) isequal(d, indices_test), {d1, d2, d3, d4, d5, d6})));

rms_axes    = sqrt(mean(resid .^ 2, 2))';  % [Fx Fy Fz Tx Ty Tz]
labels_axes = {'', '', '', 'n_x', 'n_y', 'n_z'};
units_axes  = {'', '', '', 'Nm', 'Nm', 'Nm'};
range_all   = (max(w_hand, [], 2) - min(w_hand, [], 2))';

cx_final = phi_final(2) / m_known;  cy_final = phi_final(3) / m_known;  cz_final = phi_final(4) / m_known;
cx_cad   = phi_sw(2) / phi_sw(1);   cy_cad   = phi_sw(3) / phi_sw(1);   cz_cad   = phi_sw(4) / phi_sw(1);

log_msg('run_regression', '');
log_msg('run_regression', 'Inertial - CAD prior in italics');
print_param('m',   '[kg]',      m_known,  phi_sw(1));
print_param('c_x', '[m]',       cx_final, cx_cad);
print_param('c_y', '[m]',       cy_final, cy_cad);
print_param('c_z', '[m]',       cz_final, cz_cad);
print_param('Ixx', '[kg m^2]',  phi_final(5),  phi_sw(5));
print_param('Iyy', '[kg m^2]',  phi_final(8),  phi_sw(8));
print_param('Izz', '[kg m^2]',  phi_final(10), phi_sw(10));
print_param('Ixy', '[kg m^2]',  phi_final(6),  phi_sw(6));
print_param('Ixz', '[kg m^2]',  phi_final(7),  phi_sw(7));
print_param('Iyz', '[kg m^2]',  phi_final(9),  phi_sw(9));

log_msg('run_regression', '');
log_msg('run_regression', 'Friction');
log_msg('run_regression', '%-6s %-10s %10.4f', 'bv',   '[Nm s/rad]', bv);
log_msg('run_regression', '%-6s %-10s %10.4f', 'bc',   '[Nm]',       bc);
log_msg('run_regression', '%-6s %-10s %10.4f', 'bpsi', '[Nm s/rad]', bpsi);
log_msg('run_regression', '%-6s %-10s %10.4f', 'mu_r', '[-]',        mu_r);

log_msg('run_regression', '');
log_msg('run_regression', 'Regression - *effective kappa');
log_msg('run_regression', '%-6s %-10s %10s', 'Train', '[sec]', sprintf('%.0f', training_data.N * dt));
log_msg('run_regression', '%-6s %-10s %10s', 'Test',  '[sec]', sprintf('%.0f', testing_data.N * dt));
log_msg('run_regression', '%-6s %-10s %10s', 'lambda', '[-]',  sprintf('%.0f', lambda));
log_msg('run_regression', '%-6s %-10s %10s', 'kappa*', '[-]',  sprintf('%.1f', kappa_eff));

log_msg('run_regression', '');
log_msg('run_regression', 'RMS residuals - %% is percent of signal range');
for i = 4:6
    pct = 100 * rms_axes(i) / range_all(i);
    log_msg('run_regression', '%-6s %-10s %10.4f', labels_axes{i}, sprintf('[%s]', units_axes{i}), rms_axes(i));
    log_msg('run_regression', '%-6s %-10s %10.1f', '', '[%]', pct);
end

%%

log_msg('run_regression', 'plotting test-set results');

fig = figure('Color', 'w', 'Name', 'Test Results', 'Units', 'inches', 'Position', [0 0 5 5.8]);
t_col_test = t_col_all(indices_test);

yaw_unwrapped = unwrap(data.yaw(indices_test));
yaw_all = rad2deg(yaw_unwrapped - yaw_unwrapped(1));

c_red = [228 26 28] / 255;
c_blue = [55 126 184] / 255;
c_blk = [0 0 0];
c_grn = [77 175 74] / 255;

pad = 0.35;
tq_idx = [4 5 6];
tq_labels = {'$x$ Torque (Nm)', '$y$ Torque (Nm)', '$z$ Torque (Nm)'};

tl = tiledlayout(5, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
legs = gobjects(5, 1);
axs = gobjects(5, 1);

axs(1) = nexttile;
plot(t_col_test, data.r_WB_W(1, indices_test), 'Color', c_red, 'LineWidth', 1.0); hold on;
plot(t_col_test, data.r_WB_W(2, indices_test), 'Color', c_grn, 'LineWidth', 1.0);
ylabel('Position (m)', 'Interpreter', 'latex'); grid on; grid minor;
legs(1) = legend({'$x$', '$y$'}, 'Interpreter', 'latex', 'Location', 'southeast', ...
    'FontSize', 6, 'Orientation', 'horizontal');
legs(1).BoxFace.ColorType = 'truecoloralpha';
legs(1).BoxFace.ColorData = uint8([255; 255; 255; 180]);
ylim_pad([data.r_WB_W(1, indices_test), data.r_WB_W(2, indices_test)], pad);
xlim([t_col_test(1), t_col_test(end)]);
set(gca, 'XTickLabel', [], 'FontSize', 7, 'TickLabelInterpreter', 'latex', 'Box', 'off');

axs(2) = nexttile;
plot(t_col_test, rad2deg(data.pitch(indices_test)), 'Color', c_grn, 'LineWidth', 1.0); hold on;
plot(t_col_test, rad2deg(data.roll(indices_test)), 'Color', c_red, 'LineWidth', 1.0);
plot(t_col_test, yaw_all, 'Color', c_blue, 'LineWidth', 1.0);
ylabel('Orientation (deg)', 'Interpreter', 'latex'); grid on; grid minor;
legs(2) = legend({'Pitch', 'Camber', 'Yaw'}, 'Interpreter', 'latex', 'Location', 'southeast', ...
    'FontSize', 6, 'Orientation', 'horizontal');
legs(2).BoxFace.ColorType = 'truecoloralpha';
legs(2).BoxFace.ColorData = uint8([255; 255; 255; 180]);
ylim_pad([rad2deg(data.pitch(indices_test)), rad2deg(data.roll(indices_test)), yaw_all], pad);
xlim([t_col_test(1), t_col_test(end)]);
set(gca, 'XTickLabel', [], 'FontSize', 7, 'TickLabelInterpreter', 'latex', 'Box', 'off');

for j = 1:3
    axs(j + 2) = nexttile;
    i = tq_idx(j);
    d1 = abs(w_hand(i, :) - w_bal(i, :));
    d2 = w_bal(i, :);
    d3 = w_hand(i, :);
    plot(t_col_test, d1, '-', 'Color', c_blue, 'LineWidth', 0.8); hold on;
    plot(t_col_test, d2, '-', 'Color', c_red, 'LineWidth', 1.0);
    plot(t_col_test, d3, '-', 'Color', c_blk, 'LineWidth', 1.0);
    ylabel(tq_labels{j}, 'Interpreter', 'latex');
    grid on; grid minor;
    legs(j + 2) = legend({'$|$Res$|$', 'Est', 'Ref'}, 'Interpreter', 'latex', 'Location', 'southeast', ...
        'FontSize', 6, 'Orientation', 'horizontal');
    legs(j + 2).BoxFace.ColorType = 'truecoloralpha';
    legs(j + 2).BoxFace.ColorData = uint8([255; 255; 255; 180]);
    ylim_pad([d1, d2, d3], pad);
    xlim([t_col_test(1), t_col_test(end)]);
    set(gca, 'FontSize', 7, 'TickLabelInterpreter', 'latex', 'Box', 'off');

    if j < 3
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 7);
    end

end

drawnow;
repin_all(legs, axs);
fig.SizeChangedFcn = @(~, ~) repin_all(legs, axs);

function [res, J] = evaluate_cost_function(x, pre, m_known, x_prior, lambda)
    res = compute_residuals(x, pre, m_known, x_prior, lambda);

    if nargout > 1
        J = compute_jacobian(x, pre, m_known, lambda);
    end

end

function res = compute_residuals(x, d, m_known, x_prior, lambda)

    N = d.N;

    [phi, bv, bc, b_psi] = logchol_to_real(x, m_known);

    res = zeros(3 * N, 1);

    for k = 1:N
        f = evaluate_f(d, k, phi, bv, bc, b_psi);
        
        rows = (3 * k - 2):(3 * k);
        res(rows) = d.C * d.G * (d.y(:, k) - f - d.b);
    end

    theta_m = x(1:9);
    theta_m_CAD = x_prior(1:9);
    
    damp = [lambda * (theta_m - theta_m_CAD); zeros(3, 1)];
    res = [res; damp];

end

function f = evaluate_f(d, k, phi, bv, bc, b_psi)

    wC = [0; 
          0; 
          0; 
          0; 
          0; 
         -b_psi * d.psidot(k)];
    
     wA = (-bv * d.omega_w(k) - bc * sign(d.omega_w(k))) * [zeros(3, 1); 
                                                            [0;0;1]];
    
      f = d.W_B_C(:, :, k) * d.A(:, :, k) * phi - wC - d.W_B_C(:, :, k) * wA;

end

function J = compute_jacobian(x, d, m_known, lambda)
    
    N = d.N;
    J = zeros(3 * N, 12);
        
    [phi, bv, bc, b_psi] = logchol_to_real(x, m_known);

    dphi_dd = compute_dphi_dtheta_quotient(phi, x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), m_known);

    for k = 1:N
        
        Jk_C = [d.W_B_C(:, :, k) * d.A(:, :, k) * dphi_dd, ...
                    d.W_B_C(:, :, k) * d.omega_w(k) * [zeros(3, 1); 0;0;1] * bv, ...
                    d.W_B_C(:, :, k) * sign(d.omega_w(k)) * [zeros(3, 1); 0;0;1] * bc, ...
                    [0; 0; 0; 0; 0; d.psidot(k)] * b_psi];


        rows = (3 * k - 2):(3 * k);
        J(rows, :) = -d.C * d.G * Jk_C;
    end

    J_damp = zeros(12, 12);
    J_damp(1:9, 1:9) = lambda * eye(9);
    J = [J; J_damp];
end

function [phi, bv, bc, b_psi] = logchol_to_real(x, m_known)
    
    theta_m = x(1:9);
    theta_f = x(10:12);
    
    alpha = 0.5 * log(m_known / (theta_m(7) ^ 2 + theta_m(8) ^ 2 + theta_m(9) ^ 2 + 1));
    
    phi = theta_to_phi([alpha; theta_m]);
    
    bv    = exp(theta_f(1));
    bc    = exp(theta_f(2));
    b_psi = exp(theta_f(3));
end

function [phi_final, bv, bc, bpsi] = extract_wheelbarrow_results(x_opt, m_known)

    t1 = x_opt(7);
    t2 = x_opt(8);
    t3 = x_opt(9);

    alpha = 0.5 * log(m_known / (t1 ^ 2 + t2 ^ 2 + t3 ^ 2 + 1));

    phi_final = theta_to_phi([alpha; x_opt(1:9)]);

    bv   = exp(x_opt(10));
    bc   = exp(x_opt(11));
    bpsi = exp(x_opt(12));

end

function pre = get_data_selection(data, indices)

    pre = struct();

    pre.N = length(indices);
    pre.C = data.C;
    pre.G = data.G;
    pre.b = data.b;
    pre.A = data.A(:, :, indices);
    pre.y = data.y(:, indices);
    pre.omega_w = data.omega_w(indices);
    pre.psidot = data.psidot(indices);
    pre.z_B_C = data.z_B_C(:, indices);
    pre.W_B_C = data.W_B_C(:, :, indices);
end

function dphi_dd = compute_dphi_dtheta_quotient(phi, d1, d2, d3, s12, s23, s13, t1, t2, t3, m_known)

    ed1 = exp(d1);
    ed2 = exp(d2);
    ed3 = exp(d3);
    D   = t1^2 + t2^2 + t3^2 + 1;
    mD  = m_known / D;      
    mD2 = m_known / D^2;    

    dn = zeros(10, 9);

    dn(:,1) = [   0; ed1*t1;      0;      0;       0; -ed1*s12; -ed1*s13; 2*ed1^2;        0; 2*ed1^2];
    dn(:,2) = [   0;      0; ed2*t2;      0; 2*ed2^2;        0;        0;       0; -ed2*s23; 2*ed2^2];
    dn(:,3) = [   0;      0;      0; ed3*t3; 2*ed3^2;        0;        0; 2*ed3^2;        0;       0];
    dn(:,4) = [   0;      0;     t1;      0;   2*s12;     -ed1;        0;       0;     -s13;   2*s12];
    dn(:,5) = [   0;      0;      0;     t2;   2*s23;        0;        0;   2*s23;     -ed2;       0];
    dn(:,6) = [   0;      0;      0;     t1;   2*s13;        0;     -ed1;   2*s13;     -s12;       0];
    dn(:,7) = [2*t1;    ed1;    s12;    s13;       0;        0;        0;       0;        0;       0];
    dn(:,8) = [2*t2;      0;    ed2;    s23;       0;        0;        0;       0;        0;       0];
    dn(:,9) = [2*t3;      0;      0;    ed3;       0;        0;        0;       0;        0;       0];

    dphi_dd = mD * dn;
    dphi_dd(:,7) = dphi_dd(:,7) - (2*t1/D) * phi;
    dphi_dd(:,8) = dphi_dd(:,8) - (2*t2/D) * phi;
    dphi_dd(:,9) = dphi_dd(:,9) - (2*t3/D) * phi;

end

function print_param(name, unit, value, cad_value)
    log_msg('run_regression', '%-6s %-10s %10.4f', name, unit, value);
    log_msg('run_regression', '%-6s %-10s %10.4f', '', '', cad_value);
end
