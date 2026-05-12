clear;
clc;


% 选择基底
% basis = 'Heisenberg-Weyl';
% basis = 'Weyl';
basis = 'Pauli';

% 选择要检测纠缠的态
% state = 'state_bipartite_ppt';
state = 'state_bipartite_bounded';
% state = 'state_bipartite_horodecki';

% 储存最优的参数
p_optimal = 1;                  % 检测为纠缠的概率 p
l_optimal = 1;                  % l = length(u) = length(v)
x_optimal = 1;                  % 向量 u = (x,...,x)^{T}
function_uv = @(u) u.^(1/2);    % 设置 u 和 v 的关系 v = f(u)
alpha_optimal = [];             % 向量 α
beta_optimal = [];              % 向量 β

% 根据量子态，设置遍历范围及确定维度
if strcmp(state, 'state_bipartite_ppt')
    p_range = linspace(0.8823, 0.8, 1000);
    l_range = 1:5;
    m_range = 1:4;  % m = length(alpha) = length(beta)
    x_range = linspace(1, 0.001, 500);
    d = [3, 3];
elseif strcmp(state, 'state_bipartite_bounded')
    p_range = linspace(0.3, 0, 1000);
    l_range = 2:4;
    m_range = 2:4;  % m = length(alpha) = length(beta)
    x_range = linspace(1, 0, 501);
    d = [2, 4];
elseif strcmp(state, 'state_bipartite_horodecki')
    p_range = linspace(1, 0, 501);
    l_range = 1:4;
    m_range = 1:4;  % m = length(alpha) = length(beta)
    x_range = linspace(1, 0, 501);
    d = [3, 3];
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
term_d_kappa = (d.^2 - d) ./ kappa;

% 遍历寻找可检测出纠缠的 p 的范围
for p = p_range 
    fprintf('Try: p = %g\n', p);

    % 标记当前 p 的取值是否被检测出纠缠
    flag_p = 0;

    % 获取量子态的密度矩阵
    if strcmp(state, 'state_bipartite_ppt')
        rho = state_bipartite_ppt(p);
    elseif strcmp(state, 'state_bipartite_bounded')
        a = 0.9;
        rho = state_bipartite_bounded(p, a);
    elseif strcmp(state, 'state_bipartite_horodecki')
        a = 0.9;
        rho = state_bipartite_horodecki(p, a);
    else
        error('Undefined state.');
    end

    % 计算量子态在广义Bloch表示下的相关张量和约化密度矩阵
    [T_A, T_B, T_AB] = cofficient_tensor_bipartite(rho, d, basis);
    
    % 遍历 u = (x,...,x), 长度为l
    for l = l_range
        for x = x_range
            % 获取向量 u,v
            u = repmat(x, l, 1);
            v = function_uv(u);

            % 计算向量 u,v 的Frobenius范数的平方
            Fnorm2_u = norm(u, 'fro')^2;
            Fnorm2_v = norm(v, 'fro')^2;
            
            % 重复随机生成长度为m的向量 α,β，
            for m = m_range
                % 设置随机生成次数
                for rand_cout = 1:50
                    % 随机生成向量 α,β
                    alpha = randn(m, 1);
                    beta = randn(m, 1);

                    % 计算向量 α,β 的Frobenius范数的平方
                    Fnorm2_alpha = norm(alpha, 'fro')^2;
                    Fnorm2_beta = norm(beta, 'fro')^2;
                    
                    % A|B 划分下的扩展相关张量
                    M_A_B = extended_correlation_tensor_bipartite(T_A, T_B, T_AB, u, v, alpha, beta);

                    % 计算迹范数
                    tr_A_B = trace_norm_matrix(M_A_B);
                    
                    % 计算上界
                    upperbound = sqrt( (Fnorm2_u + Fnorm2_alpha * term_d_kappa(1)) * ...
                        (Fnorm2_v + Fnorm2_beta * term_d_kappa(2)) );

                    % 计算迹范数和上界的差值 trace - upperbound
                    difference = tr_A_B - upperbound;
                    
                    % 检测到纠缠，储存相应的参数
                    if difference > 0
                        fprintf('trace - upperbound = %g\n', difference);
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