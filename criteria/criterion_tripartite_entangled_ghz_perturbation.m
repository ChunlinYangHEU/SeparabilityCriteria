clc;
clear;


% 选择基底
% basis = 'Heisenberg-Weyl';
% basis = 'Weyl';
basis = 'Pauli';

% 选择要检测纠缠的态
state = 'state_tripartite_ghz_perturbation';

% 储存最优的参数
p_optimal = 1;                  % 可检测为纠缠的概率 p

% 根据量子态设置遍历范围及确定维度
if strcmp(state, 'state_tripartite_ghz_perturbation')
    p_range = linspace(0.6, 0.3, 10001);
    epsilon = 1;
    dims = [2, 2, 2];
else
    error('Undefined state.');
end

% 根据基底确定 κ 的值
if strcmp(basis, 'Weyl') || strcmp(basis, 'Heisenberg-Weyl')
    kappa = dims;
elseif strcmp(basis, 'Pauli')
    kappa = 2 * ones(size(dims));
else
    error('Undefined basis.');
end

% 计算上界
upper_bound = prod(sqrt((dims.^2-dims)./kappa));

% 遍历寻找可检测出纠缠的 p 的范围
for p = p_range 
    fprintf('Try: p = %g\n', p);

    flag_p = 0;

    % 获取量子态的密度矩阵
    if strcmp(state, 'state_tripartite_ghz_perturbation')
        rho = state_tripartite_ghz_perturbation(p, epsilon);
    else
        error('Undefined state.');
    end

    % 计算量子态在广义Bloch表示下的相关张量
    T = correlation_tensor_multipartite(rho, dims, basis);
    
    % 计算相关张量的迹范数
    tr = trace_norm_tensor_mixed_mode(T);
    
    % 计算迹范数和上界的差值
    difference = tr - upper_bound;

    % 检测到纠缠，储存相应的参数
    if difference > 0
        p_optimal = p;
        flag_p = 1;
    end

    if flag_p == 0
        break;
    end
end

fprintf('p_optimal = %g\n', p_optimal);
