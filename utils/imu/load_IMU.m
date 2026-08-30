function [a, aa, av, g_est] = load_IMU(filename, win, deg, downsample, start_index, start_offset, n_samples)

    % from CAD
    T_IMU_LC = [ 0.0877  -0.9962   0  -0.0861;
                -0.9962  -0.0877   0   0.0125;
                      0        0  -1  -0.0290;
                      0        0   0        1];

    R_IMU_LC   = T_IMU_LC(1:3, 1:3);
    r_LCIMU_LC = T_IMU_LC(1:3, 4);

    dexes = (start_index + start_offset):(start_index + start_offset + n_samples);

    log_msg('load_IMU', 'loading %s', filename);

    IMU_data = dlmread(filename, ',', 38, 1);

    log_msg('load_IMU', 'raw rows: %d, valid (nonzero) rows: %d', size(IMU_data, 1), nnz(IMU_data(:, 1) ~= 0));

    a_IMU_IMU = extract_a_IMU(IMU_data);
    av_IMU    = extract_av_IMU(IMU_data);
    g_est_IMU = extract_g_est(IMU_data);

    aa_IMU = differentiate_and_sgolay(av_IMU, 1, 0.001, 51, 3);

    % rotate to LC frame
    a_IMU_LC = R_IMU_LC * a_IMU_IMU;
    aa_LC    = R_IMU_LC * aa_IMU;
    av_LC    = R_IMU_LC * av_IMU;
    g_est_LC = R_IMU_LC * g_est_IMU;

    % compute acceleration at load cell frame
    a_LC_LC = zeros(size(aa_LC));

    for i = 1:size(a_LC_LC, 2)
        a_LC_LC(:, i) = a_IMU_LC(:, i) + cross(aa_LC(:, i), -r_LCIMU_LC) + cross(av_LC(:, i), cross(av_LC(:, i), -r_LCIMU_LC));
    end

    a_LC_LC  = differentiate_and_sgolay(a_LC_LC,  0, 0.001, win, deg);
    aa_LC    = differentiate_and_sgolay(aa_LC,    0, 0.001, win, deg);
    av_LC    = differentiate_and_sgolay(av_LC,    0, 0.001, win, deg);
    g_est_LC = differentiate_and_sgolay(g_est_LC, 0, 0.001, win, deg);

    a     = downsample_window(a_LC_LC,  downsample, dexes);
    aa    = downsample_window(aa_LC,    downsample, dexes);
    av    = downsample_window(av_LC,    downsample, dexes);
    g_est = downsample_window(g_est_LC, downsample, dexes);

    log_msg('load_IMU', 'windowed to %d samples (%d:%d)', numel(dexes), dexes(1), dexes(end));
    log_msg('load_IMU', 'done');

end

function a_IMU_IMU = extract_a_IMU(IMU_data)

    a_IMU_IMU = IMU_data(IMU_data(:, 1) ~= 0, 1:3)' * 9.81;

end

function av_IMU = extract_av_IMU(IMU_data)

    av_IMU = IMU_data(IMU_data(:, 1) ~= 0, 4:6)';

end

function g_est_IMU = extract_g_est(IMU_data)

    g_raw_IMU   = IMU_data(IMU_data(:, 7) ~= 0, [11, 9, 7])';
    n_samples   = nnz(IMU_data(:, 1) ~= 0);
    g_est_IMU   = zeros(3, n_samples);
    idx         = 1:2:size(g_est_IMU, 2);
    max_samples = min(length(idx), size(g_raw_IMU, 2));

    log_msg('load_IMU', 'gravity estimate: %d samples available, %d used', size(g_raw_IMU, 2), max_samples);

    g_est_IMU(:, idx(1:max_samples)) = g_raw_IMU(:, 1:max_samples);
    g_est_IMU(:, 2:2:end - 1)        = (g_est_IMU(:, 1:2:end - 2) + g_est_IMU(:, 3:2:end)) / 2;

end

function out = downsample_window(in, ds, dexes)

    out = in(:, 1:ds:end);

    if dexes(end) > size(out, 2)
        error('[load_IMU] requested window %d:%d exceeds downsampled length %d', ...
            dexes(1), dexes(end), size(out, 2));
    end

    out = out(:, dexes);

end
