function [Q, R] = classical_gram_schmidt_twice(A)
% Calls CGS twice (bad idea: it's better to do MGS twice)

    % Q1 R1 ~ A
    [Q1, R1] = classical_gram_schmidt(A);
    % Q2 R2 ~ Q1
    [Q2, R2] = classical_gram_schmidt(Q1);
    
    % Hence: Q2 R2 R1 ~ A
    
    Q = Q2;
    R = R2*R1;

end
