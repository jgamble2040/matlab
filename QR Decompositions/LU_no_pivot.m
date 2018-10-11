function [L,U] = LU_no_pivot(A)
    n=size(A,1);
    
    U=A;
    
    L=eye(n);
    
    for k = 1:(n-1)
        
        for j = (k+1) : n
            
            L(j,k) = U(j,k) / U(k,k);
            
            U(j,(k+1):n) = U(j,(k+1):n) - L(j,k)*U(k,(k+1):n);
            
            U(j,k)=0;
            
        end
        
    end
    
end

