% 系数张量
function [tensor_cell] = cofficient_tensor_tripartite(rho, dims, basis)
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
    kappa1 = kappa(1); 
    kappa2 = kappa(2); 
    kappa3 = kappa(3);
    d1 = dims(1); 
    d2 = dims(2); 
    d3 = dims(3);
    factor_A = d1 / kappa1;
    factor_B = d2 / kappa2;
    factor_C = d3 / kappa3;
    factor_AB = factor_A * factor_B;
    factor_AC = factor_A * factor_C;
    factor_BC = factor_B * factor_C;
    factor_ABC = factor_AB * factor_C;

    % 预计算单位矩阵
    I1 = eye(d1);
    I2 = eye(d2);
    I3 = eye(d3);

    % 初始化张量
    nA = d1^2 - 1;
    nB = d2^2 - 1;
    nC = d3^2 - 1;
    T_A = zeros(nA, 1);
    T_B = zeros(nB, 1);
    T_C = zeros(nC, 1);
    T_AB = zeros(nA, nB);
    T_AC = zeros(nA, nC);
    T_BC = zeros(nB, nC);
    T_ABC = zeros(nA, nB, nC);

    % 计算相关张量和约化密度矩阵
    for i = 1:nA
        opA = g{1}{i+1};
        termA = kron(kron(opA, I2), I3);
        T_A(i) = factor_A * trace(rho * termA);
        
        for j = 1:nB
            opB = g{2}{j+1};
            termAB = kron(kron(opA, opB), I3);
            T_AB(i, j) = factor_AB * trace(rho * termAB);

            for k = 1:nC
                opC = g{3}{k+1};
                termABC = kron(kron(opA, opB), opC);
                T_ABC(i, j, k) = factor_ABC * trace(rho * termABC);
            end
        end
        
        for k = 1:nC
            opC = g{3}{k+1};
            termAC = kron(kron(opA, I2), opC);
            T_AC(i, k) = factor_AC * trace(rho * termAC);
        end
    end
    
    for j = 1:nB
        opB = g{2}{j+1};
        termB = kron(kron(I1, opB), I3);
        T_B(j) = factor_B * trace(rho * termB);
        
        for k = 1:nC
            opC = g{3}{k+1};
            termBC = kron(kron(I1, opB), opC);
            T_BC(j, k) = factor_BC * trace(rho * termBC);
        end
    end
    
    for k = 1:nC
        opC = g{3}{k+1};
        termC = kron(kron(I1, I2), opC);
        T_C(k) = factor_C * trace(rho * termC);
    end

    tensor_cell = {{T_A, T_B, T_C}, {T_AB, T_AC, T_BC}, {T_ABC}};
end