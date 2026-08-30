clear; 
clc; 
close all;

addpath('../utils/vicon');
addpath('../utils/load_cell');
addpath('../utils/imu');
addpath('../utils/log');
addpath('../utils/filter');
addpath('../utils/rigid_body');

L1 = 11.25 * 0.0254;
L2 = 10.75 * 0.0254;
L3 = 10.25 * 0.0254;
L4 = 13.75 * 0.0254;
downsample        = 10;
start_offset      = 900;
n_samples         = 15400;
dt                = 0.01;
N                 = 15401;
VICON_start_index = 1058;
m_h     = 1.5310;
r_com_h = [-0.0007; -0.0015; 0.0260];

%============================================================================================
log_msg('build_data', 'Load VICON data');

win_vicon = 31;
deg_vicon = 4;

marker_data_s = load_VICON_raw('data/Trajectories.csv');
marker_static_data_s = load_VICON_raw('data/Static.csv');
VICON_dexes = (VICON_start_index + start_offset):(VICON_start_index + n_samples + start_offset);
marker_data = marker_data_s(:, :, VICON_dexes);

P_ref = make_static_reference(marker_static_data_s, 1, size(marker_static_data_s, 3));
align_dexes = [1,2,3,4,5,6,7,8,9,10,13,14,16,17,19,20];
aligned_marker_data = align_markers(marker_data, P_ref, align_dexes);

vicon_report(aligned_marker_data, marker_data, align_dexes);

%============================================================================================
log_msg('build_data', 'Load force/torque sensor data');

LC_start_index    = 2082;

[w_S_S, ~, t] = load_LC('data/exp1_021.mat', ...
                    downsample, ...
                    LC_start_index, ...
                    start_offset, ...
                    n_samples);


%============================================================================================
log_msg('build_data', 'Load IMU sensor data');

IMU_start_index     = 787;
win_imu             = 101;   
deg_imu             = 4;

[aS_S, aaS_S, avS_S] = load_IMU('data/IMU_Door.csv',...
                                        win_imu, ...
                                        deg_imu, ...
                                        downsample, ...
                                        IMU_start_index, ...
                                        start_offset, ...
                                        n_samples);


%============================================================================================
log_msg('build_data', 'World frames');

[R_S_W, r_WS_W] = compute_handle_frame_in_world(aligned_marker_data, [], 1, 2, 3);
[R_B_W, r_WB_W] = compute_body_frame_in_world(aligned_marker_data);


%============================================================================================
log_msg('build_data', 'handle wrench estimates');

w_S_hat = compute_handle_wrench_estimate(R_S_W, m_h, r_com_h);


%============================================================================================
log_msg('build_data', 'static preload');

static_preload_idx = 1:500;

b_S = compute_static_preload(w_S_S, w_S_hat, static_preload_idx);


%============================================================================================
log_msg('build_data', 'VICON-based dynamics');

door_angle = zeros(1, N);
for k = 1:N
    x_hat = R_B_W(1:3, 1, k);
    door_angle(k) = atan2(x_hat(2), x_hat(1));
end

door_angle = unwrap(door_angle);
door_angle = door_angle - door_angle(1);
door_angle = -door_angle;

omega_z_vicon = differentiate_and_sgolay(door_angle, 1, dt, win_vicon, deg_vicon);
alpha_z_vicon = differentiate_and_sgolay(door_angle, 2, dt, win_vicon, deg_vicon);

% -------------------- select dynamic measurement type -------------------

omega_z = omega_z_vicon;
alpha_z = alpha_z_vicon;

% omega_z = omega_z_imu;
% alpha_z = alpha_z_imu;

%============================================================================================
log_msg('build_data', 'compute tau_ext');

R_S_B    = zeros(3,3,N);
r_BS_B   = zeros(3,N);
W_S_B    = zeros(6,6,N);
wrench_B = zeros(6,N);

for k = 1:N
    R_S_B(:,:,k)  = R_B_W(:,:,k)' * R_S_W(:,:,k);
    r_BS_B(:,k)   = R_B_W(:,:,k)' * (r_WS_W(:,k) - r_WB_W(:,k));
    W_S_B(:,:,k)  = wrench_transmission(R_S_B(:,:,k), r_BS_B(:,k));
    wrench_B(:,k) = W_S_B(:,:,k) * (w_S_S(:,k) - w_S_hat(:,k) - b_S);
end

tau_ext = wrench_B(6,:);


%============================================================================================
log_msg('build_data', 'compute fourbar linkage');

[phi, nu] = freudenstein_phi(door_angle, L1, L2, L3, L4);

phi_dot = omega_z .* nu;

theta_bc = deg2rad(69);
theta_l  = deg2rad(25);
width_bc = deg2rad(2);
width_l  = deg2rad(2);

opening   = (phi_dot > 0) & (door_angle < theta_bc);
backcheck = (phi_dot > 0) & (door_angle >= theta_bc);
sweep     = (phi_dot < 0) & (door_angle > theta_l);
latch     = (phi_dot < 0) & (door_angle <= theta_l);

w_bc = 0.5*(1 + tanh((door_angle - theta_bc) ./ width_bc));
w_l = 0.5*(1 + tanh((theta_l - door_angle) ./ width_l));
w_s = 1 - w_l;


%% Quasistatic spring (closer) model and regression reference y
%============================================================================================
log_msg('build_data', 'compute closer spring model');

idx_qs_raw  = 1:6023;
deg_qs      = 3;
vel_thresh  = 0.02;   % rad/s

is_truly_static = abs(omega_z(idx_qs_raw)) < vel_thresh;
idx_qs = idx_qs_raw(is_truly_static);

nu_qs    = nu(idx_qs);
tau_s_qs = -tau_ext(idx_qs)' ./ nu_qs';

p_phi = polyfit(phi(idx_qs)', tau_s_qs, deg_qs);

tau_s_hat = polyval(p_phi, phi(:)) .* nu(:);   % nu*tau_s_hat(phi)


log_msg('build_data', 'regression reference, y');
%============================================================================================

y = tau_ext(:) + tau_s_hat;   % y = tau_ext + nu*tau_s_hat(phi)

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

data = struct();
data.N     = N;
data.t     = t;

data.theta   = door_angle;    % rad, VICON-derived
data.omega_z = omega_z; % rad/s, chosen dynamics source
data.alpha_z = alpha_z; % rad/s^2, chosen dynamics source (VICON-derived)

data.phi     = phi;           % rad, pinion angle from four-bar kinematics
data.phi_dot = phi_dot;       % rad/s, pinion angular rate
data.nu      = nu;            % dphi/dtheta, four-bar velocity ratio

data.theta_bc = theta_bc;         % rad, backcheck zone boundary
data.theta_l  = theta_l;          % rad, latch zone boundary
data.width_bc = width_bc;         % rad, backcheck zone tanh-blend width
data.width_l  = width_l;          % rad, latch zone tanh-blend width
data.opening   = opening;
data.backcheck = backcheck;
data.sweep     = sweep;
data.latch     = latch;

data.w_bc = w_bc;
data.w_s  = w_s;
data.w_l  = w_l;

data.idx_qs_raw     = idx_qs_raw;
data.is_truly_static = is_truly_static;
data.idx_qs          = idx_qs;
data.tau_s_qs         = tau_s_qs;

data.tau_ext = tau_ext;   % Nm, about hinge z-axis
data.y         = y;           % regression reference: tau_ext - nu*tau_s_hat(phi)
data.p_phi     = p_phi;       % quasistatic cubic spring coeffs (vs phi)

save('results/door_data.mat', 'data');
log_msg('build_data', 'saved results/door_data.mat (%d samples)', N);


%%
% helper functions


function [phi, dphi_dtheta] = freudenstein_phi(theta, L1, L2, L3, L4)
    K1 = L4 / L1;
    K4 = L4 / L2;
    K5 = (L3^2 - L4^2 - L1^2 - L2^2) / (2 * L1 * L2);

    D = cos(theta) - K1 + K4 * cos(theta) + K5;
    E = -2 * sin(theta);
    F = K1 + (K4 - 1) * cos(theta) + K5;

    gamma = 2 * atan((-E - sqrt(E .^ 2 - 4 * D .* F)) ./ (2 * D));
    phi = pi + theta - gamma;

    dphi_dtheta = gradient(phi) ./ gradient(theta);
end


function [R_B_W, r_WB_W] = compute_body_frame_in_world(aligned_marker_data)
    N = size(aligned_marker_data, 3);
    R_B_W  = zeros(3,3,N);
    r_WB_W = zeros(3,N);

    for k = 1:N
        h1 = squeeze(aligned_marker_data(13, :, k))';
        h2 = squeeze(aligned_marker_data(14, :, k))';
        h3 = squeeze(aligned_marker_data(16, :, k))';
        h4 = squeeze(aligned_marker_data(17, :, k))';
        h5 = squeeze(aligned_marker_data(19, :, k))';
        h6 = squeeze(aligned_marker_data(20, :, k))';

        hinge_pts = [h1'; h2'; h3'; h4'; h5'; h6'];
        mean_pt   = mean(hinge_pts, 1);
        [~, ~, V] = svd(hinge_pts - mean_pt, 0);
        z_hat     = V(:, 1);
        if z_hat(3) > 0
            z_hat = -z_hat;
        end

        O = mean_pt';

        % marker 11 is on the other side of the door...
        lc_pt = squeeze(aligned_marker_data(11, :, k))';
        x_raw = lc_pt - O;
        x_hat = x_raw - dot(x_raw, z_hat) * z_hat;
        x_hat = x_hat / norm(x_hat);
        y_hat = cross(z_hat, x_hat);
        y_hat = y_hat / norm(y_hat);

        R_B_W(:,:,k) = [x_hat, y_hat, z_hat];
        r_WB_W(:,k)  = O;
    end
end