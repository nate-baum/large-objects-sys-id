function T_align = compute_point_cloud_alignment(P_ref, marker_data)

    n_frames = size(marker_data, 3);
    T_align  = zeros(4, 4, n_frames);

    for i = 1:n_frames

        P_meas = marker_data(:, :, i)';
        valid  = ~any(isnan(P_meas), 1);

        if sum(valid) < 3
            T_align(:, :, i) = NaN;
            continue
        end

        P_meas_valid = P_meas(:, valid);
        P_ref_valid  = P_ref(:, valid);

        centroid_meas = mean(P_meas_valid, 2);
        centroid_ref  = mean(P_ref_valid, 2);

        P_meas_centered = P_meas_valid - centroid_meas;
        P_ref_centered  = P_ref_valid - centroid_ref;

        H = P_ref_centered * P_meas_centered';
        [U, ~, V] = svd(H);
        R = V * U';

        if det(R) < 0
            % SVD can return an improper rotation (reflection); flip the
            % smallest-singular-value axis to force a proper rotation.
            V(:, 3) = -V(:, 3);
            R = V * U';
        end

        t = centroid_meas - R * centroid_ref;

        T = eye(4);
        T(1:3, 1:3) = R;
        T(1:3, 4)   = t;

        T_align(:, :, i) = T;

    end

end
