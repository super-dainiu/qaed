addpath('/home/ys792/project/qaed/qtfm')
tol = eps;

order_n = [64, 128, 256, 512, 1024, 2048, 4096, 8192];
for n = order_n
    % load the matrix
    load(sprintf('skew_aed_n%d.mat', n));
    uni_err = rel_err(eye(n), Q' * Q);
    schur_err = rel_err(D, Q' * A * Q);
    eigvec_err = vec_err(A, Q, D);
    fprintf('using aed, n = %d, uni_err = %.2e, schur_err = %.2e, eigvec_err = %.2e\n', n, uni_err, schur_err, eigvec_err);

    % load the matrix
    load(sprintf('skew_iqr_n%d.mat', n));
    uni_err = rel_err(eye(n), Q' * Q);
    schur_err = rel_err(D, Q' * A * Q);
    eigvec_err = vec_err(A, Q, D);
    fprintf('using iqr, n = %d, uni_err = %.2e, schur_err = %.2e, eigvec_err = %.2e\n', n, uni_err, schur_err, eigvec_err); 
end

function err = rel_err(A, B)
err = norm(A - B, 'fro') / max(norm(A, 'fro'), norm(B, 'fro'));
end

function err = vec_err(A, P, D)
err = norm(A * P - P * D, 'fro') / ((norm(A, 'fro') + norm(D, 'fro')) * norm(P, 'fro'));
end
