% Compare QR algorithms
clear;
clf;
clc;
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 10);
set(groot, 'DefaultAxesFontSize', 14);

%%


% Fit polynomials through n points in [-1, 1]
m = 251;
t = linspace(-1, 1, m)';

degrees = 0:40;

condition = zeros(size(degrees));
metric1 = zeros(numel(degrees), 5);     % A ~ QR
metric2 = zeros(numel(degrees), 5);     % Q'Q ~ I

for iter = 1 : numel(degrees)

    % The polynomial has degree d: build the Vandermonde matrix
    d = degrees(iter);
    A = zeros(m, d+1);
    for k = 0 : d
        A(:, d+1-k) = t.^k;
    end
    
    %
    condition(iter) = cond(A);
    
    % Matlab's QR (Householder)
    [Q, R] = qr(A, 0);
    metric1(iter, 1) = norm(Q*R - A) / norm(A);
    metric2(iter, 1) = norm(Q'*Q - eye(d+1));
    if iter == numel(degrees), figure(2); subplot(2, 3, 1); imagesc(log10(abs(Q'*Q))); colorbar; title('Matlab''s QR'); axis tight; axis equal; end
    
    % CGS
    [Q, R] = classical_gram_schmidt(A);
    metric1(iter, 2) = norm(Q*R - A) / norm(A);
    metric2(iter, 2) = norm(Q'*Q - eye(d+1));
    if iter == numel(degrees), figure(2); subplot(2, 3, 2); imagesc(log10(abs(Q'*Q))); colorbar; title('CGS'); axis tight; axis equal; end
    
    % MGS
    [Q, R] = modified_gram_schmidt(A);
    metric1(iter, 3) = norm(Q*R - A) / norm(A);
    metric2(iter, 3) = norm(Q'*Q - eye(d+1));
    if iter == numel(degrees), figure(2); subplot(2, 3, 3); imagesc(log10(abs(Q'*Q))); colorbar; title('MGS'); axis tight; axis equal; end
    
    % CGS twice
    [Q, R] = classical_gram_schmidt_twice(A);
    metric1(iter, 4) = norm(Q*R - A) / norm(A);
    metric2(iter, 4) = norm(Q'*Q - eye(d+1));
    if iter == numel(degrees), figure(2); subplot(2, 3, 5); imagesc(log10(abs(Q'*Q))); colorbar; title('CGS twice'); axis tight; axis equal; end
    
    % MGS twice
  [Q, R] = modified_gram_schmidt_twice(A);
    metric1(iter, 5) = norm(Q*R - A) / norm(A);
    metric2(iter, 5) = norm(Q'*Q - eye(d+1));
    if iter == numel(degrees), figure(2); subplot(2, 3, 6); imagesc(log10(abs(Q'*Q))); colorbar; title('MGS twice'); axis tight; axis equal; end
    
end

% Plot

figure(1);

subplot(1, 2, 1);
loglog(condition, metric1(:, 1), ...
       condition, metric1(:, 2), ...
       condition, metric1(:, 3), ...
       condition, metric1(:, 4), ...
       condition, metric1(:, 5));
title('Relative error on A: ||QR - A|| / ||A||');
xlabel('Condition number cond(A)');
   
legend('Matlab''s QR (Householder)', ...
       'Classical GS', ...
       'Modified GS', ...
       'Classical GS twice', ...
       'Modified GS twice', ...
       'Location', 'SouthEast');
set(legend, 'Box', 'off');
   

subplot(1, 2, 2);
loglog(condition, metric2(:, 1), ...
       condition, metric2(:, 2), ...
       condition, metric2(:, 3), ...
       condition, metric2(:, 4), ...
       condition, metric2(:, 5));
title('Loss of orthogonality: ||Q''Q - I||');
xlabel('Condition number cond(A)');
   
legend('Matlab''s QR (Householder)', ...
       'Classical GS', ...
       'Modified GS', ...
       'Classical GS twice', ...
       'Modified GS twice', ...
       'Location', 'NorthWest');
set(legend, 'Box', 'off');
