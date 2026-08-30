function L = vectorize_inertia(k)

    L = [k(1) k(2) k(3)   0    0    0
           0  k(1)   0  k(2) k(3)   0
           0    0  k(1)   0  k(2) k(3)];
