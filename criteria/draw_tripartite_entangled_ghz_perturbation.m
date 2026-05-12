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

Z = [];

% 根据量子态设置遍历范围及确定维度
if strcmp(state, 'state_tripartite_ghz_perturbation')
    p_range = linspace(1, 0, 1001);
    epsilon_range = linspace(0, 2, 501);
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

index = 1;

% 遍历寻找可检测出纠缠的 p 的范围
for p = p_range 
    fprintf('Try: p = %g\n', p);
    
    for epsilon = epsilon_range
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
        
        % 储存数据
        Z(index) = difference;
    
        index = index + 1;
    
        % 检测到纠缠，储存相应的参数
        if difference > 0
            p_optimal = p;
        end

    end
end

fprintf('p_optimal = %g\n', p_optimal);


%% 绘图
X = repelem(p_range, length(epsilon_range));
Y = repmat(epsilon_range, 1, length(p_range));

% 定义网格范围
[grid_X, grid_Y] = meshgrid(linspace(min(X), max(X), 100), linspace(min(Y), max(Y), 100));

% 使用griddata进行插值
grid_Z = griddata(X, Y, Z, grid_X, grid_Y, 'linear');

% 绘制三维网格图
figure;
mesh(grid_X, grid_Y, grid_Z, 'FaceColor', 'y', 'EdgeColor', 'k');

% 添加坐标轴的表示
xlabel('p', 'FontSize', 14);
ylabel('\epsilon', 'FontSize', 14);
zlabel('Trace norm - Upper bound', 'FontSize', 14);

% 添加网格线
grid on; 

% 保持当前图形，用于叠加绘制
hold on;  

% 绘制 XOY 平面
Z_plane = zeros(size(grid_Z));
surf(grid_X, grid_Y, Z_plane, 'FaceColor', 'b');

% 保存图片
if strcmp(state, 'state_tripartite_ghz_perturbation')
    view(-37.5, 60)
    print('example6', '-dpng', '-r300');
else
    error('Undefined state.');
end