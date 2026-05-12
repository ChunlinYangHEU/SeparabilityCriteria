function [rho] = state_tripartite(p)
    ket0 = [1;0;0];
    ket1 = [0;1;0];
    ket2 = [0;0;1];
    phi = (1 / sqrt(6)) * (kron(kron(ket0, ket0), ket1) + kron(kron(ket0, ket1), ket0) ...
        + kron(kron(ket1, ket0), ket0) + kron(kron(ket1, ket1), ket2) ...
        + kron(kron(ket1, ket2), ket1) + kron(kron(ket2, ket1), ket1));
    rho = ((1 - p) / 27) * eye(27) + p * (phi * phi');
end