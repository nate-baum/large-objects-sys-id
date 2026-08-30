function dphi_dtheta = compute_dphi_dtheta(theta)
% COMPUTE_DPHI_DTHETA  Jacobian of phi = theta_to_phi(theta) w.r.t. theta (10x10).
%
%   Rows:    [m, hx, hy, hz, Ixx, Ixy, Ixz, Iyy, Iyz, Izz]
%   Columns: [a, d1, d2, d3, s12, s23, s13, t1, t2, t3]

    a   = theta(1);
    d1  = theta(2);
    d2  = theta(3);
    d3  = theta(4);
    s12 = theta(5);
    s23 = theta(6);
    s13 = theta(7);
    t1  = theta(8);
    t2  = theta(9);
    t3  = theta(10);

    e2a = exp(2 * a);
    ed1 = exp(d1);
    ed2 = exp(d2);
    ed3 = exp(d3);

    phi = theta_to_phi(theta);

    dphi_dtheta = zeros(10, 10);

    % d/da: every term carries a common exp(2a) factor
    dphi_dtheta(:, 1) = 2 * phi;

    % d/dd1
    dphi_dtheta(2, 2) =  e2a * t1 * ed1;      % hx
    dphi_dtheta(6, 2) = -e2a * s12 * ed1;     % Ixy
    dphi_dtheta(7, 2) = -e2a * s13 * ed1;     % Ixz
    dphi_dtheta(8, 2) =  2 * e2a * ed1 ^ 2;   % Iyy
    dphi_dtheta(10, 2) = 2 * e2a * ed1 ^ 2;   % Izz

    % d/dd2
    dphi_dtheta(3, 3) =  e2a * t2 * ed2;      % hy
    dphi_dtheta(5, 3) =  2 * e2a * ed2 ^ 2;   % Ixx
    dphi_dtheta(9, 3) = -e2a * s23 * ed2;     % Iyz
    dphi_dtheta(10, 3) = 2 * e2a * ed2 ^ 2;   % Izz

    % d/dd3
    dphi_dtheta(4, 4) = e2a * t3 * ed3;       % hz
    dphi_dtheta(5, 4) = 2 * e2a * ed3 ^ 2;    % Ixx
    dphi_dtheta(8, 4) = 2 * e2a * ed3 ^ 2;    % Iyy

    % d/ds12
    dphi_dtheta(3, 5)  =  e2a * t1;           % hy
    dphi_dtheta(5, 5)  =  2 * e2a * s12;      % Ixx
    dphi_dtheta(6, 5)  = -e2a * ed1;          % Ixy
    dphi_dtheta(9, 5)  = -e2a * s13;          % Iyz
    dphi_dtheta(10, 5) =  2 * e2a * s12;      % Izz

    % d/ds23
    dphi_dtheta(4, 6) = e2a * t2;             % hz
    dphi_dtheta(5, 6) = 2 * e2a * s23;        % Ixx
    dphi_dtheta(8, 6) = 2 * e2a * s23;        % Iyy
    dphi_dtheta(9, 6) = -e2a * ed2;           % Iyz

    % d/ds13
    dphi_dtheta(4, 7) = e2a * t1;             % hz
    dphi_dtheta(5, 7) = 2 * e2a * s13;        % Ixx
    dphi_dtheta(7, 7) = -e2a * ed1;           % Ixz
    dphi_dtheta(8, 7) = 2 * e2a * s13;        % Iyy
    dphi_dtheta(9, 7) = -e2a * s12;           % Iyz

    % d/dt1
    dphi_dtheta(1, 8) = 2 * e2a * t1;         % m
    dphi_dtheta(2, 8) = e2a * ed1;            % hx
    dphi_dtheta(3, 8) = e2a * s12;            % hy
    dphi_dtheta(4, 8) = e2a * s13;            % hz

    % d/dt2
    dphi_dtheta(1, 9) = 2 * e2a * t2;         % m
    dphi_dtheta(3, 9) = e2a * ed2;            % hy
    dphi_dtheta(4, 9) = e2a * s23;            % hz

    % d/dt3
    dphi_dtheta(1, 10) = 2 * e2a * t3;        % m
    dphi_dtheta(4, 10) = e2a * ed3;           % hz

end
