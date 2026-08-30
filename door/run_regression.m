clear;
clc;
close all;

addpath('../utils/plot');
addpath('../utils/log');

load('results/door_data.mat', 'data');

indices_train = [6024:10833, 14129:15000];
indices_test  = 10833:14129;

training_data = get_data_selection(data, indices_train);
testing_data  = get_data_selection(data, indices_test);
log_msg('run_regression', 'training set: %d samples, testing set: %d samples', ...
    training_data.N, testing_data.N);

x0 = zeros(6,1);

opts = optimoptions('lsqnonlin', ...
    'Algorithm',                'trust-region-reflective', ...
    'MaxIterations',            1000, ...
    'MaxFunctionEvaluations',   50000, ...
    'FunctionTolerance',        1e-12, ...
    'StepTolerance',            1e-12, ...
    'Display',                  'iter');

cost = @(x) evaluate_cost_function(x, training_data);

log_msg('run_regression', 'starting optimization...');
[x_opt, resnorm, ~, ~, ~, ~, J_sol] = lsqnonlin(cost, x0, [], [], opts);
log_msg('run_regression', 'optimization finished, resnorm=%.4g', resnorm);

res_test = evaluate_cost_function(x_opt, testing_data);
rms_test = sqrt(mean(res_test .^ 2));

s_J = svd(full(J_sol));
kappa_J = s_J(1) / s_J(end);

f_hat = zeros(data.N, 1);
for k = 1:data.N
    f_hat(k) = evaluate_f(data, k, exp(x_opt(1)), exp(x_opt(2)), exp(x_opt(3)), exp(x_opt(4)), exp(x_opt(5)), exp(x_opt(6)));
end

log_msg('run_regression', 'theta_bc = %.1f deg   theta_l = %.1f deg   width_bc = %.1f deg   width_l = %.1f deg', ...
    rad2deg(data.theta_bc), rad2deg(data.theta_l), rad2deg(data.width_bc), rad2deg(data.width_l));

%% Per-zone RMS breakdown

log_msg('run_regression', '');
log_msg('run_regression', '=== PER-ZONE TEST RMS ===');
log_msg('run_regression', '%-12s %10s %14s', 'Zone', 'N pts', 'RMS (Nm)');
log_msg('run_regression', '%-12s %10d %14.4f', 'opening',   sum(data.opening(indices_test)),   sqrt(mean(res_test(data.opening(indices_test)).^2)));
log_msg('run_regression', '%-12s %10d %14.4f', 'backcheck', sum(data.backcheck(indices_test)), sqrt(mean(res_test(data.backcheck(indices_test)).^2)));
log_msg('run_regression', '%-12s %10d %14.4f', 'sweep',     sum(data.sweep(indices_test)),     sqrt(mean(res_test(data.sweep(indices_test)).^2)));
log_msg('run_regression', '%-12s %10d %14.4f', 'latch',     sum(data.latch(indices_test)),     sqrt(mean(res_test(data.latch(indices_test)).^2)));

%% Paper-style report table

dt_door = data.t(2) - data.t(1);
range_nz = max(data.y(indices_test)) - min(data.y(indices_test));
pct_nz   = 100 * rms_test / range_nz;

c0 = data.p_phi(4);  c1 = data.p_phi(3);  c2 = data.p_phi(2);  c3 = data.p_phi(1);

log_msg('run_regression', '');
log_msg('run_regression', 'Inertial');
log_msg('run_regression', '%-6s %-10s %10.4f', 'Izz', '[kg m^2]', exp(x_opt(1)));

log_msg('run_regression', '');
log_msg('run_regression', 'Damping');
log_msg('run_regression', '%-6s %-10s %10.4f', 'b_bc', '[Nm s/rad]', exp(x_opt(2)));
log_msg('run_regression', '%-6s %-10s %10.4f', 'b_s',  '[Nm s/rad]', exp(x_opt(3)));
log_msg('run_regression', '%-6s %-10s %10.4f', 'b_l',  '[Nm s/rad]', exp(x_opt(4)));

log_msg('run_regression', '');
log_msg('run_regression', 'Friction');
log_msg('run_regression', '%-6s %-10s %10.4f', 'b_v', '[Nm s/rad]', exp(x_opt(5)));
log_msg('run_regression', '%-6s %-10s %10.4f', 'b_c', '[Nm]',       exp(x_opt(6)));

log_msg('run_regression', '');
log_msg('run_regression', 'Stiffness');
log_msg('run_regression', '%-6s %-10s %10.4f', 'c0', '[Nm]',        c0);
log_msg('run_regression', '%-6s %-10s %10.4f', 'c1', '[Nm/rad]',    c1);
log_msg('run_regression', '%-6s %-10s %10.4f', 'c2', '[Nm/rad^2]',  c2);
log_msg('run_regression', '%-6s %-10s %10.4f', 'c3', '[Nm/rad^3]',  c3);

log_msg('run_regression', '');
log_msg('run_regression', 'Regression - *effective kappa');
log_msg('run_regression', '%-6s %-10s %10s', 'Train', '[sec]', sprintf('%.0f', training_data.N * dt_door));
log_msg('run_regression', '%-6s %-10s %10s', 'Test',  '[sec]', sprintf('%.0f', testing_data.N * dt_door));
log_msg('run_regression', '%-6s %-10s %10s', 'kappa*', '[-]', sprintf('%.1f', kappa_J));

log_msg('run_regression', '');
log_msg('run_regression', 'RMS residuals - %% is percent of signal range');
log_msg('run_regression', '%-6s %-10s %10.4f', 'n_z', '[Nm]', rms_test);
log_msg('run_regression', '%-6s %-10s %10.1f', '', '[%]', pct_nz);

t_test = data.t(indices_test);

tau_pred = f_hat(indices_test)';
tau_meas = data.y(indices_test)';
resid    = abs(res_test)';

%% Effective hydraulic damper damping vs door angle (0 to 80 deg, positive convention)

theta_max_deg = 80;
theta_max = deg2rad(theta_max_deg);

theta_open_half  = linspace(0, theta_max, 250)';    % 0 -> 80 (positive = open)
theta_close_half = linspace(theta_max, 0, 250)';    % 80 -> 0

bw_open = 0.5*(1 + tanh((theta_open_half - data.theta_bc) ./ data.width_bc));
b_eff_open_half = exp(x_opt(2)) * bw_open;

latch_w = 0.5*(1 + tanh((data.theta_l - theta_close_half) ./ data.width_l));
sweep_w = 1 - latch_w;
b_eff_close_half = exp(x_opt(3))*sweep_w + exp(x_opt(4))*latch_w;

figure('Name', 'Effective Damping', 'Color', 'w', ...
       'Units', 'inches', 'Position', [1 1 3.5 2]);
hold on;

h1 = plot(rad2deg(theta_open_half),  b_eff_open_half,  'b-', 'LineWidth', 1.0);
h2 = plot(rad2deg(theta_close_half), b_eff_close_half, 'r-', 'LineWidth', 1.0);
h3 = xline(rad2deg(data.theta_bc), '--', 'Color', 'b', 'LineWidth', 1.0);
h4 = xline(rad2deg(data.theta_l),  '--', 'Color', 'r', 'LineWidth', 1.0);

xlabel('$\theta$ (deg)', 'Interpreter', 'latex', 'FontSize', 8);
ylabel('Damping coeff (Nm s/rad)', 'Interpreter', 'latex', 'FontSize', 8);
legend([h1 h2 h3 h4], {'Opening', 'Closing', '$\theta_{bc}$', '$\theta_l$'}, ...
       'Interpreter', 'latex', 'FontSize', 7, ...
       'Orientation', 'horizontal', 'Location', 'southoutside');

grid on; grid minor;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 7);
xlim([0 theta_max_deg]);
ylim([0 max(ylim)]);

%% Quasistatic spring fit result plot (velocity-filtered)

phi_fit_range = linspace(min(data.phi(data.idx_qs)), max(data.phi(data.idx_qs)), 200)';
tau_fit_range = polyval(data.p_phi, phi_fit_range);

excluded_idx = data.idx_qs_raw(~data.is_truly_static);

figure('Name', 'Spring Fit Validation', 'Color', 'w', ...
       'Units', 'inches', 'Position', [1 1 3.5 2]);
hold on;

nu_excluded = data.nu(excluded_idx);
plot(rad2deg(data.phi(excluded_idx)), -data.tau_ext(excluded_idx)' ./ nu_excluded', '.', ...
     'Color', [0.7 0.7 0.7], 'MarkerSize', 4);
plot(rad2deg(data.phi(data.idx_qs)), data.tau_s_qs, 'k.', 'MarkerSize', 4);
plot(rad2deg(phi_fit_range), tau_fit_range, 'r-', 'LineWidth', 1.0);

xlabel('$\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 8);
ylabel('$\tau_{s}(\phi)$ (Nm)', 'Interpreter', 'latex', 'FontSize', 8);
legend('Excluded ($|\dot\theta|>$thresh)', 'Quasistatic data', 'Cubic fit', ...
       'Interpreter', 'latex', 'FontSize', 7, ...
       'Orientation', 'horizontal', 'Location', 'southoutside');

grid on; grid minor;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 7);
xlim([min(rad2deg(data.phi(data.idx_qs_raw))), max(rad2deg(data.phi(data.idx_qs_raw)))]);

%% Door test validation (two-tile: angle + torque)

fs_door     = 8;    % font size
lw_door     = 0.75; % line width
lw_res_door = 0.75; % residual line width

c_red  = [228 26  28] / 255;
c_blue = [55  126 184] / 255;
c_blk  = [0   0   0 ];

contact_test  = abs(tau_meas) > 0;
tau_pred_plot = tau_pred; tau_pred_plot(~contact_test) = NaN;
resid_plot    = abs(resid); resid_plot(~contact_test)  = NaN;

fig_door = figure('Name', 'Door Test Validation', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 5 2.5]);
pad_door = 0.25;
tl_door  = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
legs_d   = gobjects(2,1);
axs_d    = gobjects(2,1);

axs_d(1) = nexttile;
plot(t_test, rad2deg(data.theta(indices_test)), '-', 'Color', c_blk, 'LineWidth', lw_door);
ylabel('$\theta$ (deg)', 'Interpreter','latex', 'FontSize', fs_door);
grid on; grid minor;
xlim([t_test(1), t_test(end)]);
set(gca,'XTickLabel',[],'FontSize',fs_door,'TickLabelInterpreter','latex','Box','off');
legs_d(1) = legend({'$\theta$'}, 'Interpreter','latex', 'FontSize',fs_door, ...
    'Orientation','horizontal', 'Location','southeast');
legs_d(1).BoxFace.ColorType = 'truecoloralpha';
legs_d(1).BoxFace.ColorData = uint8([255; 255; 255; 180]);
ylim_pad(rad2deg(data.theta(indices_test)), pad_door);

axs_d(2) = nexttile;
plot(t_test, resid, '-', 'Color', c_blue, 'LineWidth', lw_res_door); hold on;
plot(t_test, tau_pred, '-', 'Color', c_red,  'LineWidth', lw_door);
plot(t_test, tau_meas, '-', 'Color', c_blk,  'LineWidth', lw_door);
ylabel('Torque (Nm)', 'Interpreter','latex', 'FontSize', fs_door);
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', fs_door);
xlim([t_test(1), t_test(end)]);
grid on; grid minor;
set(gca,'FontSize',fs_door,'TickLabelInterpreter','latex','Box','off');
legs_d(2) = legend({'$|$Res$|$','Est','Ref'}, 'Interpreter','latex', 'FontSize',fs_door, ...
    'Orientation','horizontal', 'Location','southeast');
legs_d(2).BoxFace.ColorType = 'truecoloralpha';
legs_d(2).BoxFace.ColorData = uint8([255; 255; 255; 180]);
ylim_pad([resid(:); tau_pred(:); tau_meas(:)], pad_door);

drawnow;
repin_all(legs_d, axs_d);
fig_door.SizeChangedFcn = @(~,~) repin_all(legs_d, axs_d);

%%

function res = evaluate_cost_function(x, d)
    I    = exp(x(1));
    b_bc = exp(x(2));
    b_s  = exp(x(3));
    b_l  = exp(x(4));
    b_v  = exp(x(5));
    b_c  = exp(x(6));

    res = zeros(d.N, 1);
    for k = 1:d.N
        f = evaluate_f(d, k, I, b_bc, b_s, b_l, b_v, b_c);
        res(k) = d.y(k) - f;
    end
end

function f = evaluate_f(d, k, I, b_bc, b_s, b_l, b_v, b_c)

    phi_dot_k = d.phi_dot(k);
    nu_k      = d.nu(k);

    if phi_dot_k > 0
        tau_b = d.w_bc(k) * b_bc * phi_dot_k;
    elseif phi_dot_k < 0
        tau_b = (d.w_s(k) * b_s + d.w_l(k) * b_l) * phi_dot_k;
    else
        tau_b = 0;
    end

    f = I * d.alpha_z(k) + b_v * d.omega_z(k) + b_c * sign(d.omega_z(k)) + nu_k * tau_b;
end

function pre = get_data_selection(data, indices)

    pre = struct();

    pre.N                = length(indices);
    pre.alpha_z          = data.alpha_z(indices);
    pre.omega_z          = data.omega_z(indices);
    pre.phi_dot          = data.phi_dot(indices);
    pre.nu               = data.nu(indices);
    pre.w_bc = data.w_bc(indices);
    pre.w_s  = data.w_s(indices);
    pre.w_l  = data.w_l(indices);
    pre.y                = data.y(indices);
end
