% Compute the trace norm of a matix
function [max_trace_norm] = trace_norm_tensor_mixed_mode(A)
    N = ndims(A);

    max_trace_norm = 0;

    for k = 1:N-1
        % 生成所有大小为k的R子集
        R_subsets = nchoosek(1:N, k);

        for i = 1:size(R_subsets, 1)
            R = R_subsets(i, :);
            C = setdiff(1:N, R);
            
            % 遍历所有可能的n和m
            for n = 1:k
                for m = 1:N-k
                    % 获取展开后的矩阵 A_{(R,n;C,m)}
                    matrix = matricization_mixed_mode(A, R, n, C, m);
                    
                    % 计算展开矩阵的迹范数
                    current_trace_norm = trace_norm_matrix(matrix);
                    
                    % 更新最大值
                    if current_trace_norm > max_trace_norm
                        max_trace_norm = current_trace_norm;
                    end
                end
            end
        end
    end
end