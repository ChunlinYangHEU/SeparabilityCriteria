% 生成带噪声的PPT纠缠态
function [rho] = state_bipartite_ppt(p)
    % 预计算
    identity = eye(9);
    ket0 = [1; 0; 0];
    ket1 = [0; 1; 0];
    ket2 = [0; 0; 1];
    psi0 = (1 / sqrt(2)) * kron(ket0, ket0 - ket1);
    psi1 = (1 / sqrt(2)) * kron(ket0 - ket1, ket2);
    psi2 = (1 / sqrt(2)) * kron(ket2, ket1 - ket2);
    psi3 = (1 / sqrt(2)) * kron(ket1 - ket2, ket0);
    psi4 = (1 / 3) * kron(ket0 + ket1 + ket2, ket0 + ket1 + ket2);
    
    % 构造PPT纠缠态
    rho = (1 / 4) * (identity - psi0 * (psi0') - psi1 * (psi1') - psi2 * (psi2') - psi3 * (psi3') - psi4 * (psi4'));

    % 加上噪声
    rho = ((1 - p) / 9) * identity + p * rho;
end