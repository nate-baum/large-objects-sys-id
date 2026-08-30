function ylim_pad(data_vals, pad)

    data_vals = data_vals(~isnan(data_vals) & isfinite(data_vals));
    if isempty(data_vals); return; end
    dmin = min(data_vals);
    dmax = max(data_vals);
    rng = dmax - dmin;
    if rng == 0; ylim([dmin - 1, dmax + 1]); return; end
    ylim([dmin - pad * rng, dmax + 0.05 * rng]);

end
