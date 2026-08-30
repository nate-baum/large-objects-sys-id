clear; 
clc; 
close all;

addpath('../utils/vicon');
addpath('../utils/load_cell');
addpath('../utils/imu');
addpath('../utils/log');
addpath('../utils/filter');
addpath('../utils/rigid_body');

downsample   = 10;
start_offset = 1;
n_samples    = 21200;
N            = 21201;
dt           = 0.01;
m_h          = 0.8082;
r_com_h      = [0.0002; 0.0001; 0.0814];
g_W          = [0; 0; -9.81];
R_w          = 0.197; 
win_vicon    = 11; 
deg          = 4;

%============================================================================================
log_msg('build_data', 'Load CAD prior');

phi_sw = build_sw_prior();


%============================================================================================
log_msg('build_data', 'Load VICON data');

VICON_start_index  = 1008;
static_ref_start   = 20951;
static_ref_end     = 21201;

[marker_data, marker_map] = load_VICON_raw('data/WB Trajectories.csv');
VICON_dexes = (VICON_start_index + start_offset):(VICON_start_index + start_offset + n_samples);
marker_data = marker_data(:, :, VICON_dexes);

P_ref = make_static_reference(marker_data, static_ref_start, static_ref_end);
aligned_marker_data = align_markers(marker_data, P_ref);

vicon_report(aligned_marker_data, marker_data);


%============================================================================================
log_msg('build_data', 'Load force/torque sensor data');

LC_start_index     = 2005;

[w_R_R, w_L_L, t] = load_LC('data/exp1_039.mat', ...
                        downsample, ...
                        LC_start_index, ...
                        start_offset, ... 
                        n_samples);


%============================================================================================
log_msg('build_data', 'Load IMU sensor data');

IMU_start_index    = 1439;
win_imu            = 101;

[a_R_R, aa_R_R, av_R_R] = load_IMU('data/WB_Object_IMU.csv', ...
                              win_imu, ...
                              deg, ...
                              downsample, ...
                              IMU_start_index, ...
                              start_offset, ...
                              n_samples);


%============================================================================================
log_msg('build_data', 'World frames');

[R_B_W, r_WB_W] = compute_body_frame_in_world(aligned_marker_data, marker_map);
[R_L_W, r_WL_W] = compute_handle_frame_in_world(aligned_marker_data, marker_map, 'WBObject_WB_ID14', 'WBObject_WB_ID15', 'WBObject_WB_ID16');
[R_R_W, r_WR_W] = compute_handle_frame_in_world(aligned_marker_data, marker_map, 'WBObject_WB_ID5', 'WBObject_WB_ID4', 'WBObject_WB_ID3');


%============================================================================================
log_msg('build_data', 'handle wrench estimates');

w_L_hat = compute_handle_wrench_estimate(R_L_W, m_h, r_com_h);
w_R_hat = compute_handle_wrench_estimate(R_R_W, m_h, r_com_h);


%============================================================================================
log_msg('build_data', 'static preload');

static_preload_idx = 1000:1500;

b_R = compute_static_preload(w_R_R, w_R_hat, static_preload_idx);
b_L = compute_static_preload(w_L_L, w_L_hat, static_preload_idx);


%============================================================================================
log_msg('build_data', 'IMU-based dynamics');

R_R_B      = zeros(3,3,N);
r_BR_B     = zeros(3,N);
av_B_imu   = zeros(3,N);
aa_B_imu   = zeros(3,N);
a_B_B_imu  = zeros(3,N);
for k = 1:N
    
    R_R_B(:,:,k) = R_B_W(:,:,k)' * R_R_W(:,:,k);
    r_BR_B(:,k)  = R_B_W(:,:,k)' * (r_WR_W(:,k) - r_WB_W(:,k));

    av_B_imu(:,k) = R_R_B(:,:,k) * av_R_R(:,k);
    aa_B_imu(:,k) = R_R_B(:,:,k) * aa_R_R(:,k);

    a_R_B = R_R_B(:,:,k) * a_R_R(:,k);
    
    a_B_B_imu(:,k)  = a_R_B - cross(aa_B_imu(:,k), r_BR_B(:,k)) ...
        - cross(av_B_imu(:,k), cross(av_B_imu(:,k), r_BR_B(:,k)));
end


%============================================================================================
log_msg('build_data', 'VICON-based dynamics');

a_B_W_vicon = differentiate_and_sgolay(r_WB_W, 2, dt, win_vicon, deg);

a_B_B_vicon = zeros(3,N);
for k = 1:N
    a_B_B_vicon(:,k) = R_B_W(:,:,k)' * a_B_W_vicon(:,k) - R_B_W(:,:,k)' * g_W;
end

av_B_vicon_raw = compute_angular_velocity_from_rotations(R_B_W, dt);

av_B_vicon = differentiate_and_sgolay(av_B_vicon_raw, 0, dt, win_vicon, deg);
aa_B_vicon = differentiate_and_sgolay(av_B_vicon_raw, 1, dt, win_vicon, deg);

% -------------------- select dynamic measurement type -------------------

av_B  = av_B_vicon;    % or av_B_imu
aa_B  = aa_B_vicon;    % or aa_B_imu
a_B_B = a_B_B_vicon;   % or a_B_B_imu


%============================================================================================
log_msg('build_data', 'psidot');

psidot = zeros(1,N);
for k = 1:N
    av_W = R_B_W(:,:,k) * av_B(:,k);
    psidot(k) = av_W(3);
end


%============================================================================================
log_msg('build_data', 'regressor matrix A');

A = zeros(6,10,N);
for k = 1:N
    a_ddg = a_B_B(:,k);
    Sw    = skew(av_B(:,k));
    Swdot = skew(aa_B(:,k));
    Lw    = vectorize_inertia(av_B(:,k));
    Lwdot = vectorize_inertia(aa_B(:,k));

    A(:,:,k) = [a_ddg,       Swdot+Sw*Sw,  zeros(3,6);
                zeros(3,1),  skew(-a_ddg), Lwdot+Sw*Lw];
end


%============================================================================================
log_msg('build_data', 'contact frame in world frame A');

r_WC_W = zeros(3,N);
for k = 1:N
    r_WC_W(:,k) = [r_WB_W(1,k); r_WB_W(2,k); 0];
end
R_C_W = compute_contact_to_world_rotation(R_B_W);


%============================================================================================
log_msg('build_data', 'transmission matrix from B to C');

R_B_C = zeros(3,3,N);
r_CB_C = zeros(3,N);
W_B_C = zeros(6,6,N);
for k = 1:N
    R_B_C(:,:,k) = R_C_W(:,:,k)' * R_B_W(:,:,k);
    r_CB_C(:,k) = R_C_W(:,:,k)' * (r_WB_W(:,k) - r_WC_W(:,k));
    W_B_C(:,:,k) = wrench_transmission(R_B_C(:,:,k), r_CB_C(:,k));
end


%============================================================================================
log_msg('build_data', 'transmission matrix from R to C');

W_R_C = zeros(6,6,N);
for k = 1:N
    R_R_C = R_C_W(:,:,k)' * R_R_W(:,:,k);
    r_CR_C = R_C_W(:,:,k)' * (r_WR_W(:,k) - r_WC_W(:,k));
    W_R_C(:,:,k) = wrench_transmission(R_R_C, r_CR_C);
end


%============================================================================================
log_msg('build_data', 'transmission matrix from L to C');

W_L_C = zeros(6,6,N);
for k = 1:N
    R_L_C = R_C_W(:,:,k)' * R_L_W(:,:,k);
    r_CL_C = R_C_W(:,:,k)' * (r_WL_W(:,k) - r_WC_W(:,k));
    W_L_C(:,:,k) = wrench_transmission(R_L_C, r_CL_C);
end


%============================================================================================
log_msg('build_data', 'axle axis in C frame');

z_B_C = zeros(3,N);
for k = 1:N
    z_B_C(:,k) = R_C_W(:,:,k)' * R_B_W(:,:,k) * [0;0;1];
end


%============================================================================================
log_msg('build_data', 'orientation (pitch/roll/yaw)');

pitch = zeros(1,N);
roll  = zeros(1,N);
yaw   = zeros(1,N);
for k = 1:N
    pitch(k) = asin(R_B_W(3,2,k));
    roll(k)  = asin(R_B_W(3,3,k));
    x_C_W = R_C_W(:,1,k);
    yaw(k) = atan2(x_C_W(2), x_C_W(1));
end


%============================================================================================
log_msg('build_data', 'contact velocity in W frame');

v_C_W = differentiate_and_sgolay(r_WC_W, 1, dt, win_vicon, deg);
v_C_C = zeros(3,N);
omega_w = zeros(1,N);
for k = 1:N
    v_C_C(:,k) = R_C_W(:,:,k)' * v_C_W(:,k);
    omega_w(k) = v_C_C(1,k) / R_w;
end


%============================================================================================
log_msg('build_data', 'regression reference, y');

y = zeros(6,N);
for k = 1:N
    y(:,k) = W_R_C(:,:,k) * (w_R_R(:,k) - w_R_hat(:,k) - b_R) ...
           + W_L_C(:,:,k) * (w_L_L(:,k) - w_L_hat(:,k) - b_L);
end

% -------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
data.N       = N;
data.A       = A;
data.y       = y;
data.omega_w = omega_w;
data.psidot  = psidot;
data.z_B_C   = z_B_C;
data.W_B_C   = W_B_C;
data.r_WB_W  = r_WB_W;
data.pitch = pitch;
data.roll  = roll;
data.yaw   = yaw;

save('results/wheelbarrow_data.mat', 'data', 'phi_sw', 'dt');
log_msg('build_data', 'saved wheelbarrow_data.mat (%d samples)', N);


%%
% helper functions

function [R_B_W, r_WB_W] = compute_body_frame_in_world(aligned_marker_data, marker_map)
    p_1_W  = get_marker(aligned_marker_data, marker_map, 'WBObject_WB_ID8');
    p_2_W  = get_marker(aligned_marker_data, marker_map, 'WBObject_WB_ID11');
    p_a1_W = get_marker(aligned_marker_data, marker_map, 'WBObject_WB_ID9');
    p_a2_W = get_marker(aligned_marker_data, marker_map, 'WBObject_WB_ID6');
    p_b1_W = get_marker(aligned_marker_data, marker_map, 'WBObject_WB_ID10');
    p_b2_W = get_marker(aligned_marker_data, marker_map, 'WBObject_WB_ID13');

    N = size(p_1_W, 2);
    R_B_W  = zeros(3, 3, N);
    r_WB_W = zeros(3, N);
    for k = 1:N
        z_raw = p_2_W(:,k) - p_1_W(:,k);
        y_raw = 0.5*((p_a2_W(:,k) - p_a1_W(:,k)) + (p_b2_W(:,k) - p_b1_W(:,k)));
        z_hat = z_raw / norm(z_raw);
        x_hat = cross(y_raw, z_hat);  x_hat = x_hat / norm(x_hat);
        y_hat = cross(z_hat, x_hat);
        R_B_W(:,:,k)  = [x_hat, y_hat, z_hat];
        r_WB_W(:,k)   = 0.5*(p_1_W(:,k)+p_2_W(:,k));
    end
end

function av_raw = compute_angular_velocity_from_rotations(R, dt)
    N      = size(R, 3);
    av_raw = zeros(3, N);
    for k = 2:N-1
        R_dot = (R(:,:,k+1) - R(:,:,k-1)) / (2*dt);
        S = R(:,:,k)' * R_dot;
        av_raw(:,k) = [S(3,2); S(1,3); S(2,1)];
    end
    av_raw(:,1)   = av_raw(:,2);
    av_raw(:,end) = av_raw(:,N-1);
end

function R_C_W = compute_contact_to_world_rotation(R_B_W)
    N     = size(R_B_W, 3);
    R_C_W = zeros(3, 3, N);
    z_C_W = [0; 0; 1];
    for k = 1:N
        z_B_W = R_B_W(:,3,k);
        y_C_W = (z_B_W - (z_C_W' * z_B_W) * z_C_W);
        y_C_W = y_C_W / norm(y_C_W);
        x_C_W = cross(y_C_W, z_C_W);
        R_C_W(:,:,k) = [x_C_W, y_C_W, z_C_W];
    end
end

function phi_sw = build_sw_prior()
    m     = 45.0;
    r_com = [-0.1121; 0.4498; 0.0004];
    I_com = [3.6670  0.0553  0.0046;
             0.0553  1.0924  0.0008;
             0.0046  0.0008  3.1350];
    phi_sw = com_inertia_to_phi(m, r_com, I_com);
end
