%% Setup a matrix A and a vector b for a
%  poorly-conditioned  least-squares problem

clear; close all; clc;

m = 251;
t = linspace(-1, 1, m)';

d = 21;
A = zeros(m, d+1);
for k = 0 : d
    A(:, d+1-k) = t.^k;
end

x = randn(d+1, 1);
b = A*x;

fprintf('Condition number of A: %.4e\n\n', cond(A));

%% Naive approach: solve the normal equations
x_nrm = (A'*A)\(A'*b);
err_x = norm(x - x_nrm)/norm(x);
err_b = norm(b - A*x_nrm)/norm(b);
fprintf('Normal equations\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

%% Using Matlab's backslash (which uses a QR with pivoting)
x_bck = A\b;
err_x = norm(x - x_bck)/norm(x);
err_b = norm(b - A*x_bck)/norm(b);
fprintf('Matlab''s backslash (QR with pivot)\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

% Notice how we get about cond(A) error reduction in both metrics

%% Using Classical Gram-Schmidt QR algorithm
[Q, R] = classical_gram_schmidt(A);
x_cgs = R\(Q'*b);
err_x = norm(x - x_cgs)/norm(x);
err_b = norm(b - A*x_cgs)/norm(b);
fprintf('Classical Gram-Schmidt\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

% Classical Gram-Schmidt fails miserably... But it's not because QR is a
% bad choice: see next with Matlab's QR algorithm.

%% Using Matlab's QR algorithm
[Q, R] = qr(A, 0);
x_qr = R\(Q'*b);
err_x = norm(x - x_qr)/norm(x);
err_b = norm(b - A*x_qr)/norm(b);
fprintf('Matlab''s QR algorithm\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

% Matlab's backslash uses QR with pivoting -- already with simple QR, we
% get most of the benefits.

%% Using Modified Gram-Schmidt QR algorithm
[Q, R] = modified_gram_schmidt(A);
x_mgs = R\(Q'*b);
err_x = norm(x - x_mgs)/norm(x);
err_b = norm(b - A*x_mgs)/norm(b);
fprintf('Modified Gram-Schmidt\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

% Modified Gram-Schmidt helps a lot, but it's not quite there yet.

%% Using Modified Gram-Schmidt twice QR algorithm
[Q, R] = modified_gram_schmidt_twice(A);
x_mgs2 = R\(Q'*b);
err_x = norm(x - x_mgs2)/norm(x);
err_b = norm(b - A*x_mgs2)/norm(b);
fprintf('Modified Gram-Schmidt twice\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

% MGS twice is as good as Matlab's QR (but is more expensive
% computationally.)

%% Using MGS on [A, b]
[QQ, RR] = modified_gram_schmidt([A, b]);
r = RR(1:end-1, end);
R = RR(1:end-1, 1:end-1);
x_qr2 = R\r;
err_x = norm(x - x_qr2)/norm(x);
err_b = norm(b - A*x_qr2)/norm(b);
fprintf('Modified Gram-Schmidt on [A, b]\n\tRelative error on x: %.4e\n\tRelative residue:    %.4e\n\n', err_x, err_b);

% This formulation relies less on orthonormality of Q, hence performance is
% less degraded by the lack of orthonormality due to numerical round-off.
