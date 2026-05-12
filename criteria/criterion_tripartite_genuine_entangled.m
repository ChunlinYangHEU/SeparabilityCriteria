clc;
clear;


% 选择基底
basis = 'Heisenberg-Weyl';
% basis = 'Weyl';
% basis = 'Pauli';

% 选择要检测纠缠的态
state = 'state_tripartite_werner_type';
% state = 'state_tripartite';

% 储存最优的参数
p_optimal = 1;                  % 可检测为纠缠的概率 p
l_optimal = 1;                  % l = length(u) = length(v)
x_optimal = 1;                  % 向量 u = (x,...,x)^{T}
function_uv = @(u) u.^(1/2);    % 设置 u 和 v 的关系 v = f(u)
alpha_optimal = [];             % 向量 α
beta_optimal = [];              % 向量 β

% 根据量子态设置遍历范围及确定维度
if strcmp(state, 'state_tripartite_werner_type')
    p_range = linspace(0.3, 0.2, 1000);
    l_range = 2:4;
    m_range = 2:4;  % m = length(alpha) = length(beta)
    x_range = linspace(1, 0.001, 1000);
    d = [3, 3, 2];
elseif strcmp(state, 'state_tripartite')
    p_range = linspace(1, 0, 501);
    l_range = 1:4;
    m_range = 1:4;  % m = length(alpha) = length(beta)
    x_range = linspace(1, 0.001, 500);
    d = [3, 3, 3];
else
    error('Undefined state.');
end

% 根据基底确定 κ 的值
if strcmp(basis, 'Weyl') || strcmp(basis, 'Heisenberg-Weyl')
    kappa = d;
elseif strcmp(basis, 'Pauli')
    kappa = 2 * ones(size(d));
else
    error('Undefined basis.');
end

% 预计算
m_A = upper_bound_correlation_tensor(d(1), kappa(1));
m_B = upper_bound_correlation_tensor(d(2), kappa(2));
m_C = upper_bound_correlation_tensor(d(3), kappa(3));
m_AB = upper_bound_correlation_tensor([d(1),d(2)], [kappa(1),kappa(2)]);
m_BC = upper_bound_correlation_tensor([d(2),d(3)], [kappa(2),kappa(3)]);
m_AC = upper_bound_correlation_tensor([d(1),d(3)], [kappa(1),kappa(3)]);

% 遍历寻找可检测出纠缠的 p 的范围
for p = p_range 
    fprintf('Try: p = %g\n', p);

    % 标记当前 p 的取值是否被检测出纠缠
    flag_p = 0;

    % 获取量子态的密度矩阵
    if strcmp(state, 'state_tripartite_werner_type')
        rho = state_tripartite_werner_type(p);
    elseif strcmp(state, 'state_tripartite')
        rho = state_tripartite(p);
    else
        error('Undefined state.');
    end

    % 计算量子态在广义Bloch表示下的相关张量和约化密度矩阵
    tensor_cell = cofficient_tensor_tripartite(rho, d, basis);
    T_A = tensor_cell{1}{1};
    T_B = tensor_cell{1}{2};
    T_C = tensor_cell{1}{3};
    T_AB = tensor_cell{2}{1};
    T_AC = tensor_cell{2}{2};
    T_BC = tensor_cell{2}{3};
    T_ABC = tensor_cell{3}{1};
    
    % 计算相关张量和约化密度矩阵的向量化后的维度
    dims_T_A = numel(T_A);
    dims_T_B = numel(T_B);
    dims_T_C = numel(T_C);
    dims_T_AB = numel(T_AB);
    dims_T_AC = numel(T_AC);
    dims_T_BC = numel(T_BC);

    % 遍历 u = (x,...,x), 长度为l
    for l = l_range
        for x = x_range
            % 获取向量 u,v
            u = repmat(x, l, 1);
            v = function_uv(u);

            % 计算向量 u,v 的Frobenius范数的平方
            Fnorm2_u = norm(u, 'fro')^2;
            Fnorm2_v = norm(v, 'fro')^2;
            
            % 重复随机生成向量 α,β，长度为m
            for m = m_range
                for rand_cout = 1:20
                    % 随机生成向量 α,β
                    alpha = randn(m, 1);
                    beta = randn(m, 1);

                    % 计算向量 α,β 的Frobenius范数的平方
                    Fnorm2_alpha = norm(alpha, 'fro')^2;
                    Fnorm2_beta = norm(beta, 'fro')^2;
                    
                    % 设置展开的模式
                    kR = 1;
                    kC = 1;
                    
                    % A|BC 划分下的扩展相关张量、迹范数和上界
                    M_A_BC = extended_correlation_tensor_tripartite( ...
                        T_A, T_BC, T_ABC, 1, kR, [2,3], kC, u, v, alpha, beta);
                    tr_A_BC = trace_norm_matrix(M_A_BC);
                    upperbound_A_BC = sqrt((Fnorm2_u + Fnorm2_alpha * m_A) * (Fnorm2_v + Fnorm2_beta * m_BC));

                    % AB|C 划分下的扩展相关张量、迹范数和上界
                    M_AB_C = extended_correlation_tensor_tripartite( ...
                        T_AB, T_C, T_ABC, [1,2], kR, 3, kC, u, v, alpha, beta);
                    tr_AB_C = trace_norm_matrix(M_AB_C);
                    upperbound_AB_C = sqrt((Fnorm2_u + Fnorm2_alpha * m_AB) * (Fnorm2_v + Fnorm2_beta * m_C));

                    % AC|B 划分下的扩展相关张量、迹范数和上界
                    M_AC_B = extended_correlation_tensor_tripartite( ...
                        T_AC, T_B, T_ABC, [1,3], kR, 2, kC, u, v, alpha, beta);
                    tr_AC_B = trace_norm_matrix(M_AC_B);
                    upperbound_AC_B = sqrt((Fnorm2_u + Fnorm2_alpha * m_AC) * (Fnorm2_v + Fnorm2_beta * m_B));

                    % 计算 M 和M0
                    M = (tr_A_BC + tr_AC_B + tr_AB_C) / 3;
                    M0 = max([upperbound_A_BC, upperbound_AB_C, upperbound_AC_B]);

                    % 计算迹范数和上界的差值 M - M0
                    difference = M - M0;

                    % 检测到纠缠，储存相应的参数
                    if difference > 0
                        p_optimal = p;
                        l_optimal = l;
                        x_optimal = x;
                        alpha_optimal = alpha;
                        beta_optimal = beta;
                        fprintf('p = %g, l = %g, x = %g, alpha = [', p, l, x);
                        fprintf('%g ', alpha); fprintf('], beta = [');
                        fprintf('%g ', beta); fprintf(']\n');
                        disp('----------------------------------------------------------');

                        flag_p = 1;
                        break;
                    end
                end
                
                if flag_p == 1
                    break;
                end
            end

            if flag_p == 1
                break;
            end
        end

        if flag_p == 1
            break;
        end
    end

    % 如果对于当前的 p , 检测不出纠缠，则结束 p 的循环 
    if flag_p == 0
        disp('Not found');
        disp('----------------------------------------------------------');
        break;
    end
end

fprintf('p_optimal = %g, l_optimal = %g, x_optimal = %g, alpha_optimal = [', p_optimal, l_optimal, x_optimal);
fprintf('%g ', alpha_optimal); fprintf('], beta_optimal = [');
fprintf('%g ', beta_optimal); fprintf(']\n');