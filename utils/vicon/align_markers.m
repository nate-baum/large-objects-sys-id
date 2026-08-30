function [aligned_marker_data, T_align] = align_markers(marker_data, P_ref, align_dexes)

    if nargin < 3 || isempty(align_dexes)
        align_dexes = 1:size(marker_data, 1);
    end

    n_missing = nnz(any(any(isnan(marker_data), 2), 1));

    log_msg('align_markers', '%d / %d frames have at least one missing marker', n_missing, size(marker_data, 3));

    T_align = compute_point_cloud_alignment(P_ref(:, align_dexes), marker_data(align_dexes, :, :));
    aligned_marker_data = apply_alignment(marker_data, P_ref, T_align);

end

function marker_data = apply_alignment(marker_data, P_ref, T_align)

    M = size(P_ref, 2);

    for frame_idx = 1:size(T_align, 3)

        T_current = T_align(:, :, frame_idx);

        P_projected = T_current * [P_ref; ones(1, M)];

        marker_data(:, :, frame_idx) = P_projected(1:3, :)';

    end

end
