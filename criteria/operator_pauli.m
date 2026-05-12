% Pauli算子基底
function [g] = operator_pauli(d)
    n = log2(d);
    pauli = cell(1, 4);
    pauli{1} = [1, 0;
                0, 1];
    pauli{2} = [0, 1;
                1, 0];
    pauli{3} = [0, -1i;
                1i, 0];
    pauli{4} = [1, 0;
                0, -1];
    g = cell(1, d^2);
    for i = 1:d^2
        i_bin = dec2base(i-1, 4, n);
        g{i} = 1;
        for j = i_bin
            g{i} = kron(g{i}, pauli{str2num(j)+1});
        end
    end
end