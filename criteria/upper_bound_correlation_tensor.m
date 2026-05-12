% 计算相关张量Fronenius范数平方的上界
function [m] = upper_bound_correlation_tensor(d, kappa)
    % 预计算
    N = numel(d);
    nkappa = numel(kappa);
    d = max(d);
    D = prod(d);
    K = prod(kappa);

    if N ~= nkappa
        error('The length of d and kappa should be the same.');
    end
    
    if N == 1
        m = (d^2-d) / kappa;
    elseif N == 2 || ( N>=3 && D/(d^2)<1 )
        m = (D/K) * (D + 1/(N-1) - (D/(N-1)) * sum(1./(D.^2)));
    elseif N >= 3 && D/(d^2)>=1
        m = (D/K) * (D + 2/(N-2) - (D/N-2) * sum(1./(D.^2)));
    end
end