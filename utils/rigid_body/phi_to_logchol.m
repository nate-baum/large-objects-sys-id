function theta = phi_to_logchol(phi)

    m = phi(1);
    h = phi(2:4);
    Ixx = phi(5);
    Ixy = phi(6);
    Ixz = phi(7);
    Iyy = phi(8);
    Iyz = phi(9);
    Izz = phi(10);

    I_bar = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
    Sigma = 0.5 * trace(I_bar) * eye(3) - I_bar;
    J = [Sigma, h; h', m];

    U_raw = chol(J);
    scale = U_raw(4, 4);
    alpha = log(scale);
    U_norm = U_raw / scale;

    d1 = log(U_norm(1, 1));
    d2 = log(U_norm(2, 2));
    d3 = log(U_norm(3, 3));
    s12 = U_norm(1, 2);
    s23 = U_norm(2, 3);
    s13 = U_norm(1, 3);
    t1 = U_norm(1, 4);
    t2 = U_norm(2, 4);
    t3 = U_norm(3, 4);

    theta = [alpha; d1; d2; d3; s12; s23; s13; t1; t2; t3];

end
