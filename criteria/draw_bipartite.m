clear;
clc;


% 选择要检测纠缠的态
% state = 'state_bipartite_ppt';
state = 'state_bipartite_bounded';
% state = 'state_bipartite_horodecki';

% 选择基底
% basis = 'Heisenberg-Weyl';
% basis = 'Weyl';
basis = 'Pauli';

% 设置参数
if strcmp(state, 'state_bipartite_ppt')
    function_uv = @(u) u.^(1/2);
    l = 4;
    alpha = 0.063782;
    beta = 0.0786454;
elseif strcmp(state, 'state_bipartite_bounded')
    function_uv = @(u) u.^(1/2);
    l = 2;
    alpha = [-3.23405; 1.35293];
    beta = [-1.83346; -0.969888];
elseif strcmp(state, 'state_bipartite_horodecki')
    function_uv = @(u) u.^(1/2);
    l = 2;
    alpha = [-0.797806; -0.508044; 0.0337091];
    beta = [0.800669; 0.337797; -0.362977];
else
    error('Undefined state.');
end

% 根据量子态设置遍历范围及确定维度
if strcmp(state, 'state_bipartite_ppt')
    p_start = 1;
    p_end = 0.8;
    p_step = 1000;
    x_start = 1;
    x_end = 0.001;
    x_step = 500;
    p_range = linspace(p_start, p_end, p_step);
    x_range = linspace(x_start, x_end, x_step);
    dims = [3, 3];
elseif strcmp(state, 'state_bipartite_bounded')
    p_start = 0.3;
    p_end = 0;
    p_step = 500;
    x_start = 1;
    x_end = 0;
    x_step = 500;
    p_range = linspace(p_start, p_end, p_step);
    x_range = linspace(x_start, x_end, x_step);
    dims = [2, 4];
elseif strcmp(state, 'state_bipartite_horodecki')
    p_start = 1;
    p_end = 0;
    p_step = 501;
    x_start = 1;
    x_end = 0;
    x_step = 501;
    p_range = linspace(p_start, p_end, p_step);
    x_range = linspace(x_start, x_end, x_step);
    dims = [3, 3];
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

% 计算向量 α,β 的Frobenius范数的平方
Fnorm2_alpha = norm(alpha, 'fro')^2;
Fnorm2_beta = norm(beta, 'fro')^2;

% 预计算
term_d_kappa = (dims.^2 - dims) ./ kappa;

% 存储散点数据
X = zeros(1, p_step * x_step);
Y = zeros(1, p_step * x_step);
Z = zeros(1, p_step * x_step);

index = 1;

% 遍历 p
for p = p_range 
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
    [T_A, T_B, T_AB] = cofficient_tensor_bipartite(rho, dims, basis);
    
    % 遍历 x
    for x = x_range
        % 获取向量 u,v
        u = repmat(x, l, 1);
        v = function_uv(u);

        % 计算向量 u,v 的Frobenius范数的平方
        Fnorm2_u = norm(u, 'fro')^2;
        Fnorm2_v = norm(v, 'fro')^2;
                               
        % A|B 划分下的扩展相关张量
        M_A_B = extended_correlation_tensor_bipartite(T_A, T_B, T_AB, u, v, alpha, beta);

        % 计算迹范数
        tr_A_B = trace_norm_matrix(M_A_B);
        
        % 计算上界
        upperbound = sqrt( (Fnorm2_u + Fnorm2_alpha * term_d_kappa(1)) * ...
                        (Fnorm2_v + Fnorm2_beta * term_d_kappa(2)) );

        % 计算迹范数和上界的差值 trace - upperbound
        difference = tr_A_B - upperbound;

        % 存储数据
        X(index) = x;
        Y(index) = p;
        Z(index) = difference;

        index = index + 1;
    end
end

% 定义网格范围
[grid_X, grid_Y] = meshgrid(linspace(min(X), max(X), 100), linspace(min(Y), max(Y), 100));

% 使用griddata进行插值
grid_Z = griddata(X, Y, Z, grid_X, grid_Y, 'linear');

% 绘制三维网格图
figure;
mesh(grid_X, grid_Y, grid_Z, 'FaceColor', 'y', 'EdgeColor', 'k');

% 添加坐标轴的表示
xlabel('x', 'FontSize', 14);
ylabel('p', 'FontSize', 14);
zlabel('Trace norm - Upper bound', 'FontSize', 14);

% 添加网格线
grid on; 

% 保持当前图形，用于叠加绘制
hold on;  

% 绘制 XOY 平面
Z_plane = zeros(size(grid_Z));
surf(grid_X, grid_Y, Z_plane, 'FaceColor', 'b');

% 保存图片
if strcmp(state, 'state_bipartite_ppt')
    view(-37.5+180, 60)
    print('example1', '-dpng', '-r300');
elseif strcmp(state, 'state_bipartite_bounded')
    view(-37.5, 30)
    print('example2', '-dpng', '-r300');
else
    error('Undefined state.');
end