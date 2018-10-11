function [L, U] = my_LU_no_pivot_no_force_0(A)

    n = size(A, 1);
    
    U = A;
    L = eye(n);
    
    for k = 1 : (n-1)

        for j = (k+1) : n
        
            L(j, k) = U(j, k) / U(k, k);
            
            U(j, k:n) = U(j, k:n) - L(j, k) * U(k, k:n);
            
        end

    end
    
end
