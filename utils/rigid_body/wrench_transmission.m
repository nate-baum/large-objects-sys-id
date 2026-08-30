function W = wrench_transmission(aRb, a_r_ab)

    W = [               aRb, zeros(3, 3);
         skew(a_r_ab) * aRb,        aRb];

end
