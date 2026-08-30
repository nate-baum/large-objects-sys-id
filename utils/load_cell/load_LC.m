function [w_R_R, w_L_L, t] = load_LC(filename, downsample, start_index, start_offset, n_samples)

    dexes = (start_index + start_offset):(start_index + start_offset + n_samples);

    log_msg('load_LC', 'loading %s', filename);

    data = load(filename);
    exp_name = fieldnames(data);
    X = data.(exp_name{1}).X;
    Y = data.(exp_name{1}).Y;
    [~, nY] = size(Y);

    for i = 1:nY
        path = Y(i).Path;
        path_parts = strsplit(path, '/');
        Y(i).var_names = path_parts{2};
    end

    % time vectors
    t1 = X(1).Data;
    t2 = X(3).Data;

    % get wrench from platform 1
    Y1 = Y(strcmp({Y.Device}, 'Platform') == 1);
    fx_1 = Y1(strcmp({Y1.var_names}, 'Fx1_RawZero') == 1).Data';
    fy_1 = Y1(strcmp({Y1.var_names}, 'Fy1_RawZero') == 1).Data';
    fz_1 = Y1(strcmp({Y1.var_names}, 'Fz1_RawZero') == 1).Data';
    nx_1 = Y1(strcmp({Y1.var_names}, 'Tx1_RawZero') == 1).Data';
    ny_1 = Y1(strcmp({Y1.var_names}, 'Ty1_RawZero') == 1).Data';
    nz_1 = Y1(strcmp({Y1.var_names}, 'Tz1_RawZero') == 1).Data';
    w1 = [fx_1, fy_1, fz_1, nx_1, ny_1, nz_1]';

    % get wrench from platform 2
    Y2 = Y(strcmp({Y.Device}, 'Platform_2') == 1);
    fx_2 = Y2(strcmp({Y2.var_names}, 'Fx2_RawZero') == 1).Data';
    fy_2 = Y2(strcmp({Y2.var_names}, 'Fy2_RawZero') == 1).Data';
    fz_2 = Y2(strcmp({Y2.var_names}, 'Fz2_RawZero') == 1).Data';
    nx_2 = Y2(strcmp({Y2.var_names}, 'Tx2_RawZero') == 1).Data';
    ny_2 = Y2(strcmp({Y2.var_names}, 'Ty2_RawZero') == 1).Data';
    nz_2 = Y2(strcmp({Y2.var_names}, 'Tz2_RawZero') == 1).Data';
    w2 = [fx_2, fy_2, fz_2, nx_2, ny_2, nz_2]';


    log_msg('load_LC', 'raw samples: Platform=%d, Platform_2=%d', numel(t1), numel(t2));

    % take shortest data
    minLen = min(numel(t1), numel(t2));
    t  = t1(1:minLen);
    w1 = w1(:, 1:minLen);
    w2 = w2(:, 1:minLen);


    log_msg('load_LC', 'truncated to %d aligned samples', minLen);

    t_ds  = t(1:downsample:end);
    w1_ds = w1(:, 1:downsample:end);
    w2_ds = w2(:, 1:downsample:end);


    log_msg('load_LC', 'downsampled by %d to %d samples', downsample, numel(t_ds));

    if dexes(end) > numel(t_ds)
        error('[load_LC] requested window %d:%d exceeds downsampled length %d', ...
            dexes(1), dexes(end), numel(t_ds));
    end

    t     = t_ds(dexes);
    w_R_R = w1_ds(:, dexes);
    w_L_L = w2_ds(:, dexes);

    log_msg('load_LC', 'windowed to %d samples (%d:%d)', numel(dexes), dexes(1), dexes(end));
    log_msg('load_LC', 'done');
end
