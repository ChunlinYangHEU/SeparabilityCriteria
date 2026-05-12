function [M] = matricization_mixed_mode(A, R, n, C, m)
% A   : d1×...×dN 张量
% R   : 行维度索引向量（全局编号）
% n   : 1..|R|   指 R 内部的第 n 个维度
% C   : 列维度索引向量（全局编号）
% m   : 1..|C|   指 C 内部的第 m 个维度
    
    % ---------- 检查 ----------
    assert(isequal(sort([R C]), 1:ndims(A)), ...
           'R 与 C 必须互斥且并集为 [N]');
    assert(n >=1 && n <= numel(R), 'n 超出范围');
    assert(m >=1 && m <= numel(C), 'm 超出范围');

    k = length(R);
    N = ndims(A);

    Rord = circshift(R, k-n);
    Cord = circshift(C, N-k-m);

    A = permute(A, [Rord,Cord]);

    sizeA = size(A);
    dimsR = prod(sizeA(1:k));
    dimsC = prod(sizeA(k+1:N));

    M = reshape(A, dimsR, dimsC);
end