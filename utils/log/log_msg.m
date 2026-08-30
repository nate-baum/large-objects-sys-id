function log_msg(name, fmt, varargin)

    col_width = 30;
    prefix = sprintf('[%s]', name);
    fprintf('%-*s%s\n', col_width, prefix, sprintf(fmt, varargin{:}));

end
