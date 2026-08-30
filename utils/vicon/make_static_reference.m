function P_ref = make_static_reference(marker_data, static_start, static_end)
    
    static_data = marker_data(:, :, static_start:static_end);
    
    valid_frames      = ~any(any(isnan(static_data), 2), 1);
    static_data_valid = static_data(:, :, valid_frames);
    
    P_ref = mean(static_data_valid, 3)';

    log_msg('make_static_reference', 'static reference: %d / %d valid frames', nnz(valid_frames), numel(valid_frames));

end
