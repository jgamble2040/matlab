function [Q, R] = modified_gram_schmidt_twice(A)
% Calls MGS twice

    % Q1 R1 ~ A
    [Q1, R1] = modified_gram_schmidt(A);
    % Q2 R2 ~ Q1
    [Q2, R2] = modified_gram_schmidt(Q1);
    
    % Hence: Q2 R2 R1 ~ A
    
    Q = Q2;
    R = R2*R1;

end
