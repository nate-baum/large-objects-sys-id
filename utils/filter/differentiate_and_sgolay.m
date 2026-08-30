function dx = differentiate_and_sgolay(x, p, dt, win, deg, name)

    if nargin < 6 || isempty(name)
        name = inputname(1);
    end

    if p > 0
        log_msg('differentiate_and_sgolay', '%s: differentiating (order = %d), filtering with win = %d, deg = %d', name, p, win, deg);
    else
        log_msg('differentiate_and_sgolay', '%s: filtering with win = %d, deg = %d', name, win, deg);
    end

    [~, g] = sgolay(deg, win);
    kernel = factorial(p) / (-dt) ^ p * g(:, p + 1);

    dx = zeros(size(x));

    for i = 1:size(x, 1)
        dx(i, :) = conv(x(i, :), kernel, 'same');
    end

    dx(:, 1:floor(win / 2)) = 0;
    dx(:, end - floor(win / 2):end) = 0;
