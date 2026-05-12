% 带噪声的 2×4 的有界纠缠态
function [rho] = state_bipartite_bounded(p, a)
    epsilon = (1 / sqrt(2)) * (kron([1; 0], [1; 0; 0; 0]) + kron([0; 1], [0; 1; 0; 0]));
    rho = (1 / (1 + 7 * a)) * [a, 0, 0, 0, 0, a, 0, 0;
                               0, a, 0, 0, 0, 0, a, 0;
                               0, 0, a, 0, 0, 0, 0, a;
                               0, 0, 0, a, 0, 0, 0, 0;
                               0, 0, 0, 0,       (1+a)/2, 0, 0, sqrt(1-a^2)/2;
                               a, 0, 0, 0,             0, a, 0,             0;
                               0, a, 0, 0,             0, 0, a,             0;
                               0, 0, a, 0, sqrt(1-a^2)/2, 0, 0,       (1+a)/2];
    rho = p * (epsilon * (epsilon')) + (1 - p) * rho; 
end