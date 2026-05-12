% 构造扩展相关张量
function [M] = extended_correlation_tensor_tripartite(T1, T2, T12, R, n, C, m, u, v, alpha, beta)
    
    % 预计算
    vec_T_1 = vectorization_k_mode(T1, n);
    vec_T_2 = vectorization_k_mode(T2, m);
    mat_T12 = matricization_mixed_mode(T12, R, n, C, m);
    
    % 计算每个子块
    M11 = u * (v');
    M12 = u * (kron(beta, vec_T_2)');
    M21 = kron(alpha, vec_T_1) * (v');
    M22 = kron(alpha*(beta'), mat_T12);
    
    % 将子块拼成扩展相关张量
    M = [M11, M12;
         M21, M22];
end