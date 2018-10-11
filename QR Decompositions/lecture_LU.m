% These comments are unimportant: they only tell Matlab that
% we intended to omit the semicolons (;) at the end of
% instructions. (Otherwise, the code editor issues warnings.)
%#ok<*NOPTS>
%#ok<*MNEFF>

%% On a gentle example, to check LU code is fine

clc;

A = [4 2 6 ; 2 4 1 ; 2 2 2]

[L, U] = my_LU_no_pivot(A)

A - L*U

%% On a nice matrix A, that nonetheless
%  kills the "no-pivot" algorithm
clc;

A = [0 1 ; 1 1]

[L, U] = my_LU_no_pivot(A)

A - L*U

%% Again, with "not quite 0" in A(1, 1)
clc;

A = [1e-20 1 ; 1 1]

[L, U] = my_LU_no_pivot(A)

A - L*U

%% With pivoting: everything works out.
clc;

[L, U, P] = lu(A)

P*A - L*U
