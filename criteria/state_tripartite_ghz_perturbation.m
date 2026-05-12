function rho = state_tripartite_ghz_perturbation(p, epsilon)
    I = eye(8);
    phi = (I(:,1) + epsilon * I(:,7) + I(:,8)) / sqrt(2 + epsilon^2);

    rho = ((1-p)/(8)) * I + p * (phi * (phi'));
end