% CROSSCHECK_CPP Cross-validate the C++ port against qtfm/MATLAB.
%
% Usage:
%   1. MATLAB:  generate a test matrix and save it
%        A = randq(n) .* randn(n); [~, H] = hessq(A);
%        save_qmat('H.qmat', H);
%   2. shell:   ./cpp/qaed_bench --alg aed --load H.qmat --save-out out
%   3. MATLAB:  run this script to verify out_{Q,T}.qmat against H.qmat.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'qtfm'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

H = load_qmat('H.qmat');
Q = load_qmat('out_Q.qmat');
T = load_qmat('out_T.qmat');

fprintf('||Q''Q - I||_F / sqrt(n) = %.3e\n', ...
    norm(Q' * Q - eyeq(size(Q, 1)), 'fro') / sqrt(size(Q, 1)));
fprintf('||Q''HQ - T||_F / ||H||_F = %.3e\n', ...
    norm(Q' * H * Q - T, 'fro') / norm(H, 'fro'));
fprintf('||tril(T,-1)||_F / ||T||_F = %.3e\n', ...
    norm(tril(T, -1), 'fro') / norm(T, 'fro'));

% Compare eigenvalue multisets (standardized complex eigenvalues).
d_cpp = diag(T);  ev_cpp = sort(complex(d_cpp.w + 1i * d_cpp.x));
[~, Tm] = iqrq(H); d_m = diag(Tm);
ev_m = sort(complex(d_m.w + 1i * d_m.x));
fprintf('max eigenvalue mismatch vs MATLAB iqrq = %.3e\n', ...
    max(abs(ev_cpp - ev_m)) / max(abs(ev_m)));
