% 系数张量（N体系统，N<=6）
function [T] = correlation_tensor_multipartite(rho, dims, basis)
    % 根据 basis 确定基底和 kappa 的值
    g = cell(size(dims));
    if strcmp(basis, 'Heisenberg-Weyl')
        for i1 = 1:length(dims)
            di = dims(i1);
            g{i1} = operator_hw(di);
        end
        kappa = dims;
    elseif strcmp(basis, 'Weyl')
        for i1 = 1:length(dims)
            di = dims(i1);
            g{i1} = operator_weyl(di);
        end
        kappa = dims;
    elseif strcmp(basis, 'Pauli')
        for i1 = 1:length(dims)
            di = dims(i1);
            g{i1} = operator_pauli(di);
        end
        kappa = 2 * ones(size(dims));
    else
        error('Undefined basis.');
    end

    nqubit = length(dims);

    % 预计算常数因子
    factor = prod(dims./kappa);

    % 初始化张量
    T = zeros(dims);

    % 计算相关张量
    for i1 = 1:dims(1)
        op1 = g{1}{i1+1};

        if nqubit > 1
            for i2 = 1:dims(2)
                op2 = g{2}{i2+1};

                if nqubit > 2
                    for i3 = 1:dims(3)
                        op3 = g{3}{i3+1};

                        if nqubit > 3
                            for i4 = 1:dims(4)
                                op4 = g{4}{i4+1};

                                if nqubit > 4
                                    for i5 = 1:dims(5)
                                        op5 = g{5}{i5+1};

                                        if nqubit > 5
                                            for i6 = 1:dims(6)
                                                op6 = g{6}{i6+1};
                                                op = kron(kron(kron(kron(kron(op1, op2), op3), op4), op5), op6);
                                                T(i1, i2, i3, i4, i5, i6) = factor * trace(rho * op);
                                            end
                                        else
                                            op = kron(kron(kron(kron(op1, op2), op3), op4), op5);
                                            T(i1, i2, i3, i4, i5) = factor * trace(rho * op);
                                        end
                                    end                                    
                                else
                                    op = kron(kron(kron(op1, op2), op3), op4);
                                    T(i1, i2, i3, i4) = factor * trace(rho * op);
                                end

                            end
                        else
                            op = kron(kron(op1, op2), op3);
                            T(i1, i2, i3) = factor * trace(rho * op);
                        end
                    end
                else
                    op = kron(op1, op2);
                    T(i1, i2) = factor * trace(rho * op);
                end
            end
        else
            T(i1) = factor * trace(rho * op1);
        end
    end

end

