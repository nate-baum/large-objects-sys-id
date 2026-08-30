function b = compute_static_preload(w, w_hat, idx)
    
    b = mean(w(:, idx) - w_hat(:, idx), 2);

end
