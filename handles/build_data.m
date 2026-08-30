clear;
clc;
close all;

addpath('.');
addpath('../utils/load_cell');
addpath('../utils/imu');
addpath('../utils/rigid_body');
addpath('../utils/filter');
addpath('../utils/log');

%============================================================================================
% Wheelbarrow handle

cfg = struct();
cfg.LC_filename     = 'data/exp1_036.mat';
cfg.IMU_filename    = 'data/IMU_Single.csv';
cfg.downsample      = 1;
cfg.LC_start_index  = 17757;
cfg.IMU_start_index = 9463;
cfg.start_offset    = 5000;
cfg.n_samples       = 96500;
cfg.win             = 101;
cfg.deg             = 4;
cfg.out_file        = 'results/wb_handle_data.mat';

m_h_cad = 0.9;
L_h     = 0.15;
r_h     = 0.02;

I_axial  = m_h_cad / 2  * r_h ^ 2;
I_radial = m_h_cad / 12 * (3 * r_h ^ 2 + L_h ^ 2);

cfg.m_h_cad   = m_h_cad;
cfg.I_com_cad = diag([I_radial, I_radial, I_axial]);   % cylinder axis along z
cfg.r_com_cad = [0; 0; L_h / 2];

build_handle_data(cfg);

%============================================================================================
% Door handle

cfg = struct();
cfg.LC_filename     = 'data/exp1_038.mat';
cfg.IMU_filename    = 'data/IMU_Door.csv';
cfg.downsample      = 1;
cfg.LC_start_index  = 16309;
cfg.IMU_start_index = 9760;
cfg.start_offset    = 5000;
cfg.n_samples       = 95000;
cfg.win             = 101;
cfg.deg             = 3;
cfg.out_file        = 'results/door_handle_data.mat';

m_h_cad = 1.5;
L_h     = 0.3;
r_h     = 0.02;

I_axial  = m_h_cad / 2  * r_h ^ 2;
I_radial = m_h_cad / 12 * (3 * r_h ^ 2 + L_h ^ 2);
I_com_local = diag([I_axial, I_radial, I_radial]);   % cylinder axis along x

ang = deg2rad(-45);
Rz  = [cos(ang) -sin(ang) 0;
       sin(ang)  cos(ang) 0;
       0         0        1];

cfg.m_h_cad   = m_h_cad;
cfg.I_com_cad = Rz * I_com_local * Rz';
cfg.r_com_cad = Rz * [0; 0; 0.02];

build_handle_data(cfg);

function build_handle_data(cfg)

    log_msg('build_data', 'Load force/torque sensor data');

    [w, ~, t] = load_LC(cfg.LC_filename, cfg.downsample, cfg.LC_start_index, cfg.start_offset, cfg.n_samples);


    log_msg('build_data', 'Load IMU sensor data');

    [a, aa, av, ~] = load_IMU(cfg.IMU_filename, cfg.win, cfg.deg, cfg.downsample, ...
        cfg.IMU_start_index, cfg.start_offset, cfg.n_samples);

        
    log_msg('build_data', 'tare force/torque channels');

    sigma_f = [950 * 0.015; 950 * 0.015; 1900 * 0.010];   % Fx, Fy, Fz
    sigma_n = [ 40 * 0.0125;  40 * 0.0125;  40 * 0.0125]; % Tx, Ty, Tz
    sigma   = [sigma_f; sigma_n];

    G = diag(sigma_n(1) ./ sigma);

    log_msg('build_data', 'CAD cylinder prior (m=%.2fkg)', cfg.m_h_cad);

    phi_cad = com_inertia_to_phi(cfg.m_h_cad, cfg.r_com_cad, cfg.I_com_cad);

    data    = struct();
    data.N  = numel(t);
    data.t  = t;
    data.a  = a;
    data.aa = aa;
    data.av = av;
    data.w  = w;

    data.G = G;
    data.dt         = 0.001;

    save(cfg.out_file, 'data', 'phi_cad');

    log_msg('build_data', 'saved %s (%d samples)', cfg.out_file, data.N);

end
