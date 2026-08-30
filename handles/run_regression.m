clear;
clc;
close all;

addpath('.');
addpath('../utils/rigid_body');
addpath('../utils/log');

%============================================================================================
% Wheelbarrow handle

cfg = struct();
cfg.data_file   = 'results/wb_handle_data.mat';
cfg.lambda      = 10.0;
cfg.dexes_train = 20387:46965;
cfg.dexes_test  = 51066:66693;

run_handle_regression(cfg);

%============================================================================================
% Door handle

cfg = struct();
cfg.data_file   = 'results/door_handle_data.mat';
cfg.lambda      = 20.0;
cfg.dexes_train = 2736:47446;
cfg.dexes_test  = 47446:66936;

run_handle_regression(cfg);

function run_handle_regression(cfg)

    load(cfg.data_file, 'data', 'phi_cad');
    lambda = cfg.lambda;
    dt     = data.dt;

    training_data = get_data_selection(data, cfg.dexes_train);
    testing_data  = get_data_selection(data, cfg.dexes_test);

    log_msg('run_regression', 'CAD cylinder prior (m=%.2fkg)', phi_cad(1));

    theta_cad = phi_to_logchol(phi_cad);
    x_prior = [theta_cad; zeros(6, 1)];  % zero bias prior

    theta0 = x_prior;  % start from CAD prior instead of zeros

    cost = @(x) evaluate_cost_function(x, training_data, x_prior, lambda);

    % -------------------------------------------------------------------------
    % Solve
    % -------------------------------------------------------------------------
    opts = optimoptions('lsqnonlin', ...
        'Algorithm',                'trust-region-reflective', ...
        'SpecifyObjectiveGradient', true, ...
        'MaxFunctionEvaluations',   5000, ...
        'MaxIterations',            400, ...
        'Display',                  'iter');

    log_msg('run_regression', 'starting optimization...');

    [theta_opt, resnorm] = lsqnonlin(cost, theta0, [], [], opts);
    
    
    log_msg('run_regression', 'optimization finished, resnorm=%.4g', resnorm);

    phi_est   = theta_to_phi(theta_opt(1:10));
    bias_est  = theta_opt(11:16);


    log_msg('run_regression', 'identified bias wrench:');
    log_msg('run_regression', '  f_bias: [%.4f, %.4f, %.4f] N', bias_est(1), bias_est(2), bias_est(3));
    log_msg('run_regression', '  n_bias: [%.4f, %.4f, %.4f] Nm', bias_est(4), bias_est(5), bias_est(6));

    % -------------------------------------------------------------------------
    % Conditioning
    % -------------------------------------------------------------------------
    J_data = compute_jacobian(theta_opt, training_data, 0.0);
    [~, S_data, ~] = svd(J_data, 'econ');
    s_data = diag(S_data);

    kappa_eff = sqrt(s_data(1) ^ 2 + lambda ^ 2) / sqrt(s_data(end) ^ 2 + lambda ^ 2);

    % -------------------------------------------------------------------------
    % Extract and report
    % -------------------------------------------------------------------------

    % RMS on test set
    w_pred   = testing_data.A * phi_est + repmat(bias_est, testing_data.N, 1);
    res_test = testing_data.w - w_pred;
    err_test = reshape(res_test, 6, [])';
    rms_each = sqrt(mean(err_test .^ 2));

    rms_axes    = rms_each;  % [Fx Fy Fz Tx Ty Tz]
    labels_axes = {'f_x', 'f_y', 'f_z', 'n_x', 'n_y', 'n_z'};
    units_axes  = {'N', 'N', 'N', 'Nm', 'Nm', 'Nm'};

    % Signal range per axis (max - min over test set), for RMS-as-% below
    f1_range  = max(testing_data.w_meas(1:3,:)', [], 1) - min(testing_data.w_meas(1:3,:)', [], 1);
    n1_range  = max(testing_data.w_meas(4:6,:)', [], 1) - min(testing_data.w_meas(4:6,:)', [], 1);
    range_all = [f1_range, n1_range];

    % -------------------------------------------------------------------------
    % Single-column parameter table
    
    cx_est = phi_est(2) / phi_est(1);  cy_est = phi_est(3) / phi_est(1);  cz_est = phi_est(4) / phi_est(1);
    cx_cad = phi_cad(2) / phi_cad(1);  cy_cad = phi_cad(3) / phi_cad(1);  cz_cad = phi_cad(4) / phi_cad(1);

    log_msg('run_regression', '');
    log_msg('run_regression', 'Inertial - CAD prior in italics');
    print_param('m',   '[kg]',      phi_est(1),  phi_cad(1));
    print_param('c_x', '[m]',       cx_est,      cx_cad);
    print_param('c_y', '[m]',       cy_est,      cy_cad);
    print_param('c_z', '[m]',       cz_est,      cz_cad);
    print_param('Ixx', '[kg m^2]',  phi_est(5),  phi_cad(5));
    print_param('Iyy', '[kg m^2]',  phi_est(8),  phi_cad(8));
    print_param('Izz', '[kg m^2]',  phi_est(10), phi_cad(10));
    print_param('Ixy', '[kg m^2]',  phi_est(6),  phi_cad(6));
    print_param('Ixz', '[kg m^2]',  phi_est(7),  phi_cad(7));
    print_param('Iyz', '[kg m^2]',  phi_est(9),  phi_cad(9));

    log_msg('run_regression', '');
    log_msg('run_regression', 'Regression - *effective kappa');
    log_msg('run_regression', '%-6s %-10s %10s', 'Train', '[sec]', sprintf('%.0f', training_data.N * dt));
    log_msg('run_regression', '%-6s %-10s %10s', 'Test',  '[sec]', sprintf('%.0f', testing_data.N * dt));
    log_msg('run_regression', '%-6s %-10s %10s', 'lambda', '[-]',  sprintf('%.0f', lambda));
    log_msg('run_regression', '%-6s %-10s %10s', 'kappa*', '[-]',  sprintf('%.1f', kappa_eff));

    log_msg('run_regression', '');
    log_msg('run_regression', 'RMS residuals - %% is percent of signal range');
    for i = 1:6
        pct = 100 * rms_axes(i) / range_all(i);
        log_msg('run_regression', '%-6s %-10s %10.4f', labels_axes{i}, sprintf('[%s]', units_axes{i}), rms_axes(i));
        log_msg('run_regression', '%-6s %-10s %10.1f', '', '[%]', pct);
    end

    % -------------------------------------------------------------------------
    % Plot test fit
    % -------------------------------------------------------------------------
    W_pred  = reshape(w_pred, 6, [])';
    f1_pred = -W_pred(:, 1:3);
    n1_pred = -W_pred(:, 4:6);

    f1_meas = testing_data.w_meas(1:3,:)';
    n1_meas = testing_data.w_meas(4:6,:)';

    figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 3.5 3]);
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    t_col = (0:size(f1_meas, 1) - 1) * dt;
    t_mid = t_col(end) / 2;

    nexttile;
    plot(t_col, vecnorm(f1_pred, 2, 2), 'r', 'LineWidth', 0.8); hold on;
    plot(t_col, vecnorm(f1_meas, 2, 2), 'k', 'LineWidth', 0.8);
    plot(t_col, vecnorm(f1_pred - f1_meas, 2, 2), 'b', 'LineWidth', 0.8);
    ylabel('$\|\mathbf{f}\|$ (N)', 'Interpreter', 'latex', 'FontSize', 8);
    grid on; grid minor;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 7, 'XTickLabel', []);
    xlim([t_mid - 5, t_mid + 5]);

    nexttile;
    plot(t_col, vecnorm(n1_pred, 2, 2), 'r', 'LineWidth', 0.8); hold on;
    plot(t_col, vecnorm(n1_meas, 2, 2), 'k', 'LineWidth', 0.8);
    plot(t_col, vecnorm(n1_pred - n1_meas, 2, 2), 'b', 'LineWidth', 0.8);
    ylabel('$\|\mathbf{n}\|$ (Nm)', 'Interpreter', 'latex', 'FontSize', 8);
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 8);
    grid on; grid minor;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 7);
    xlim([t_mid - 5, t_mid + 5]);
    legend('Measured', 'Predicted', 'Residual', ...
        'Interpreter', 'latex', 'FontSize', 7, ...
        'Orientation', 'horizontal', 'Location', 'southoutside');

end

function pre = get_data_selection(data, indices)

    pre = struct();
    pre.N = numel(indices);

    aa =  data.aa(:, indices);
    av =  data.av(:, indices);
    a  =  data.a(:, indices);
    f  = -data.w(1:3, indices);
    n  = -data.w(4:6, indices);

    A = zeros(6 * pre.N, 10);
    w = zeros(6 * pre.N, 1);

    for k = 1:pre.N
    
        idx   = (k-1)*6 + (1:6);
    
        Swdot = skew(aa(:,k));
        Sw    = skew(av(:,k));
        Lwdot = vectorize_inertia(aa(:,k));
        Lw    = vectorize_inertia(av(:,k));
        acc   = a(:,k);

        A(idx, :) = [         acc,  Swdot + Sw*Sw,      zeros(3,6);
                      zeros(3,1),      skew(-acc),  Lwdot + Sw*Lw];

        w(idx, 1) = [f(:,k); n(:,k)];
    
    end

    pre.A = A;
    pre.w = w;
    pre.G = repmat(diag(data.G), pre.N, 1);

    pre.t      = data.t(indices);
    pre.w_meas = data.w(:, indices);

end

function [res, J] = evaluate_cost_function(x, d, x_prior, lambda)
    res = compute_residuals(x, d, x_prior, lambda);

    if nargout > 1
        J = compute_jacobian(x, d, lambda);
    end

end

function res = compute_residuals(x, d, x_prior, lambda)

    theta  = x(1:10);
    b_bias = x(11:16);
    phi    = theta_to_phi(theta);

    r_data = d.G .* (d.w - d.A * phi - repmat(b_bias, d.N, 1));

    % DLS toward CAD prior (inertial only, not bias)
    damp = zeros(16, 1);
    damp(1:10) = lambda * (x(1:10) - x_prior(1:10));

    res = [r_data; damp];

end

function J = compute_jacobian(x, d, lambda)

    theta = x(1:10);

    J_inertial = d.A * compute_dphi_dtheta(theta);
    J_bias     = repmat(eye(6), d.N, 1);
    J_data     = -d.G .* [J_inertial, J_bias];

    J_damp = zeros(16, 16);
    J_damp(1:10, 1:10) = lambda * eye(10);

    J = [J_data; J_damp];

end

function print_param(name, unit, value, cad_value)
    log_msg('run_regression', '%-6s %-10s %10.4f', name, unit, value);
    log_msg('run_regression', '%-6s %-10s %10.4f', '', '', cad_value);
end

