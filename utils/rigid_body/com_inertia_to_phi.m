function phi = com_inertia_to_phi(m, r_com, I_com)

    h      = m * r_com;
    I_axle = I_com + m * (dot(r_com, r_com) * eye(3) - r_com * r_com');

    phi = [m; h;
           I_axle(1,1); I_axle(1,2); I_axle(1,3);
           I_axle(2,2); I_axle(2,3); I_axle(3,3)];

end
