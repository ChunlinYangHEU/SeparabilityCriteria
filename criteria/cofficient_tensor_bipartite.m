% 系数张量
function [T_A, T_B, T_AB] = cofficient_tensor_bipartite(rho, dims, basis)
    % 根据 basis 确定基底和 kappa 的值
    g = cell(size(dims));
    N = length(dims);
    if strcmp(basis, 'Heisenberg-Weyl')
        for i = 1:N
            g{i} = operator_hw(dims(i));
        end
        kappa = dims;
    elseif strcmp(basis, 'Weyl')
        for i = 1:N
            g{i} = operator_weyl(dims(i));
        end
        kappa = dims;
    elseif strcmp(basis, 'Pauli')
        for i = 1:N
            g{i} = operator_pauli(dims(i));
        end
        kappa = 2 * ones(size(dims));
    else
        error('Undefined basis.');
    end

    % 预计算常数因子
    kappaA = kappa(1); 
    kappaB = kappa(2);
    dA = dims(1); 
    dB = dims(2);
    factor_A = dA / kappaA;
    factor_B = dB / kappaB;
    factor_AB = factor_A * factor_B;

    % 预计算单位矩阵
    IA = eye(dA);
    IB = eye(dB);

    % 初始化张量
    nA = dA^2 - 1;
    nB = dB^2 - 1;
    T_A = zeros(nA, 1);
    T_B = zeros(nB, 1);
    T_AB = zeros(nA, nB);

    % 计算相关张量和约化密度矩阵
    for i = 1:nA
        opA = g{1}{i+1};
        termA = kron(opA, IB);
        T_A(i) = factor_A * trace(rho * termA);
        
        for j = 1:nB
            opB = g{2}{j+1};
            termAB = kron(opA, opB);
            T_AB(i, j) = factor_AB * trace(rho * termAB);
        end
    end
    
    for j = 1:nB
        opB = g{2}{j+1};
        termB = kron(IA, opB);
        T_B(j) = factor_B * trace(rho * termB);
    end
end