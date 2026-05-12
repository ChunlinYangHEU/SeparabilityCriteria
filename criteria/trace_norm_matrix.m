% Compute the trace norm of a matix
function [norm] = trace_norm_matrix(A)
    norm = sum(svd(A));
end