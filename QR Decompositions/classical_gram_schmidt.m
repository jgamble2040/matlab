function [Q, R] = classical_gram_schmidt(A)
% Classical Gram-Schmidt orthonormalization algorithm (unstable)

    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    
    for j = 1 : n
        
        v = A(:, j);
        
        R(1:j-1, j) = Q(:, 1:j-1)' * v;
        
        v = v - Q(:, 1:j-1) * R(1:j-1, j);
        
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
        
    end

end
