% 构造两体量子态的扩展相关张量
function [M] = extended_correlation_tensor_bipartite(T_A, T_B, T_AB, u, v, alpha, beta)  
    % 计算每个子块
    M11 = u * (v');
    M12 = u * (kron(beta, T_B)');
    M21 = kron(alpha, T_A) * (v');
    M22 = kron(alpha * (beta'), T_AB);
    
    % 将子块拼成扩展相关张量
    M = [M11, M12;
         M21, M22];
end