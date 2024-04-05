addpath('/home/ys792/qaed/qtfm')
tol = eps;

order_n = [64, 128, 256, 512, 1024, 2048, 4096, 8192];
for n = order_n
    % load the matrix
    load(sprintf('full_aed_n%d.mat', n));
    [P, D] = eigvec(Q, T);
    err = rel_err(D, inv(P) * A * P);
    fprintf('using aed, n = %d, err = %.2e\n', n, err);

    % load the matrix
    load(sprintf('full_iqr_n%d.mat', n));
    [P, D] = eigvec(Q, T);
    err = rel_err(D, inv(P) * A * P);
    fprintf('using iqr, n = %d, err = %.2e\n', n, err);
end

function err = rel_err(A, B)
err = norm(A - B, 'fro') / max(norm(A, 'fro'), norm(B, 'fro'));
end
