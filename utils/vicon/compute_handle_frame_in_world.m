function [R_W, r_WW] = compute_handle_frame_in_world(aligned_marker_data, marker_map, id1, id2, id3)

    p1 = get_marker(aligned_marker_data, marker_map, id1);
    p2 = get_marker(aligned_marker_data, marker_map, id2);
    p3 = get_marker(aligned_marker_data, marker_map, id3);

    % from CAD
    T_marker = [0.0997 -0.9950      0 -0.0897;
                0.9950  0.0997      0 -0.0201;
                     0       0 1.0000 -0.0066;
                     0       0      0 1.0000];

    N    = size(p1, 2);
    R_W  = zeros(3, 3, N);
    r_WW = zeros(3, N);

    for k = 1:N
    
        T_W_k = compute_load_cell_marker_frame(p1(:, k), p2(:, k), p3(:, k)) / T_marker;
        
        R_W(:, :, k) = T_W_k(1:3, 1:3);
        r_WW(:, k)   = T_W_k(1:3, 4);
    
    end

end

function T = compute_load_cell_marker_frame(p1, p2, p3)
    
    O1_L = p1;
    O1_M = p2;
    O1_R = p3;

    x = (O1_R - O1_M) / norm(O1_R - O1_M);
    v = (O1_L - O1_M) / norm(O1_L - O1_M);

    z = cross(v, x) / norm(cross(v, x));

    y = -cross(x, z);

    T = [x, y, z, O1_M; 0 0 0 1];

end
