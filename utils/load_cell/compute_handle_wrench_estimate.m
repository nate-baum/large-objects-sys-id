function w_hat = compute_handle_wrench_estimate(R_W, m_h, r_com_h)

    N = size(R_W, 3);
    g_W = [0; 0; -9.81];
    w_hat = zeros(6, N);

    W_grav = wrench_transmission(eye(3), r_com_h);

    for k = 1:N
        g_local = R_W(:, :, k)' * g_W;
        w_hat(:, k) = W_grav * [m_h * g_local; zeros(3, 1)];
    end

end
