%% Practice code re: solving linear systems & least squares
%% LU Factorization
% A=LU no pivot

A = [4 2 6 ; 2 4 1 ; 2 2 2];

[L, U] = LU_no_pivot(A);

A - L*U;

clc;

A = [0 1; 1 1];

[L, U] = LU_no_pivot(A);

A-L*U;
% It broke.
% with pivot

A = [0 1; 1 1];

[L, U, P] = lu(A);
%Seems like the Matlab code with pivoting.

P*A-L*U;
%Now it works
%% 
