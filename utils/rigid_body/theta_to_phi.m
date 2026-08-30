function phi = theta_to_phi(theta)

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

    m = e2a * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2 + 1);
    
    h_x = e2a * (t1 * exp(d1));
    h_y = e2a * (t1 * s12 + t2 * exp(d2));
    h_z = e2a * (t1 * s13 + t2 * s23 + t3 * exp(d3));

    Ixx = e2a * (s12 ^ 2 + exp(2 * d2) + s13 ^ 2 + s23 ^ 2 + exp(2 * d3));
    Iyy = e2a * (exp(2 * d1) + s13 ^ 2 + s23 ^ 2 + exp(2 * d3));
    Izz = e2a * (exp(2 * d1) + s12 ^ 2 + exp(2 * d2));
    Ixy = -e2a * (s12 * exp(d1));
    Iyz = -e2a * (s12 * s13 + s23 * exp(d2));
    Ixz = -e2a * (s13 * exp(d1));

    phi = [m; h_x; h_y; h_z; Ixx; Ixy; Ixz; Iyy; Iyz; Izz];

end
