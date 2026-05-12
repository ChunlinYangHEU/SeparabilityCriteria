%% 3体量子态
function [rho] = state_tripartite_werner_type(p)
    phi = (1/sqrt(5)) * ( kron(kron([0;1;0], [1;0;0]) + kron([0;0;1], [0;1;0]), [1;0]) ...
        + kron(kron([1;0;0], [1;0;0]) + kron([0;1;0], [0;1;0]) + kron([0;0;1], [0;0;1]), [0;1]) );
    rho = ((1-p)/18) * eye(18) + p * (phi * phi');
end