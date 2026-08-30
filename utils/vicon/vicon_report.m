function marker_errs = vicon_report(aligned_marker_data, marker_data, marker_indices)

    if nargin < 3 || isempty(marker_indices)
        marker_indices = 1:size(marker_data, 1);
    end

    log_msg('vicon_report', '=== Marker Reconstruction Error (mm) ===');
    log_msg('vicon_report', '%-10s %12s', 'Marker', 'Mean Err');
    log_msg('vicon_report', '%s', repmat('-', 1, 24));

    marker_errs = zeros(1, numel(marker_indices));

    for i = 1:numel(marker_indices)
        
        idx = marker_indices(i);
        
        p_est  = squeeze(aligned_marker_data(idx, :, :));
        p_meas = squeeze(marker_data(idx, :, :));
        valid  = ~any(isnan(p_meas), 1);
        
        marker_errs(i) = mean(vecnorm(p_est(:, valid) - p_meas(:, valid), 2, 1)) * 1000;

        log_msg('vicon_report', '%-10d %12.2f', idx, marker_errs(i));
    
    end

    log_msg('vicon_report', '%s', repmat('-', 1, 24));
    log_msg('vicon_report', '%-10s %12.2f', 'Mean', mean(marker_errs));

end
