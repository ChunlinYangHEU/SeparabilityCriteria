% k-mode vectorization of tensors
% vec_k(A) = vec(A_k)
function [vec_k_A] = vectorization_k_mode(A, k)
    A_k = tenmat(A, k).data;
    vec_k_A = A_k(:);
end