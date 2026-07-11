% TEST_CPP_VS_MATLAB Cross-validate and benchmark the C++ core (qaed_mex)
% against the reference MATLAB implementation.
%
% Expects qaed_mex to be built and on the path (see
% scripts/build_mex_and_check.sh). Note aedq/iqrq/hessq/eigvec now dispatch
% to the C++ core by default; qaed_accel(false) forces the pure-MATLAB
% reference paths for the baseline timings below.

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'qtfm'), fullfile(root, 'cpp'));

assert(exist('qaed_mex', 'file') == 3, 'qaed_mex not built / not on path');
rtol = eps;

fprintf('%-6s %-10s %12s %12s %10s %12s %12s\n', ...
    'n', 'alg', 't_matlab', 't_cpp', 'speedup', 'resid_cpp', 'eig_diff');

% --- small/medium n: run both, compare eigenvalues and residuals ---
for n = [64, 128, 256]
    rng(0);
    A = randq(n) .* randn(n);
    [~, H] = hessq_cpp(A);

    for alg = ["iqr", "aed"]
        qaed_accel(false);   % pure-MATLAB baseline
        if alg == "iqr"
            tic; [Qm, Tm] = iqrq(H, rtol); tm = toc;
        else
            tic; [Qm, Tm] = aedq(H, rtol); tm = toc;
        end
        qaed_accel(true);    % C++ core
        if alg == "iqr"
            tic; [Qc, Tc] = iqrq(H, rtol); tc = toc;
        else
            tic; [Qc, Tc] = aedq(H, rtol); tc = toc;
        end

        resid = norm(Qc' * H * Qc - Tc, 'fro') / norm(H, 'fro');
        unit  = norm(Qc' * Qc - eyeq(n), 'fro') / sqrt(n);
        dm = diag(Tm); dc = diag(Tc);
        em = sort(complex(dm.w + 1i * abs(dm.x)));
        ec = sort(complex(dc.w + 1i * abs(dc.x)));
        ediff = max(abs(em - ec)) / max(abs(em));

        fprintf('%-6d %-10s %12.3f %12.3f %9.1fx %12.3e %12.3e\n', ...
            n, alg + '', tm, tc, tm / tc, max(resid, unit), ediff);
    end
end

% --- larger n: C++ only (MATLAB reference would take too long) ---
qaed_accel(true);
for n = [512, 1024, 2048]
    rng(0);
    A = randq(n) .* randn(n);
    tic; [~, H] = hessq(A); th = toc;
    tic; [Qc, Tc] = aedq(H, rtol); tc = toc;
    resid = norm(Qc' * H * Qc - Tc, 'fro') / norm(H, 'fro');
    unit  = norm(Qc' * Qc - eyeq(n), 'fro') / sqrt(n);
    fprintf('%-6d %-10s %12s %12.3f %10s %12.3e (hessq %.3f s)\n', ...
        n, 'aed(cpp)', '-', tc, '-', max(resid, unit), th);
end

% --- full pipeline: patched front-end eigq and monkeypatched qtfm eig ---
n = 256; rng(0);
A = randq(n) .* randn(n);
tic; [P, D] = eigq(A, rtol); t = toc;
fprintf('eigq     n=%d: %8.3f s, ||A*P - P*D||/(||A||*||P||) = %.3e\n', ...
    n, t, norm(A * P - P * D, 'fro') / (norm(A, 'fro') * norm(P, 'fro')));

tic; [V, D2] = eig(A); t2 = toc;   % qtfm @quaternion/eig, general case
fprintf('eig(A)   n=%d: %8.3f s, ||A*V - V*D||/(||A||*||V||) = %.3e\n', ...
    n, t2, norm(A * V - V * D2, 'fro') / (norm(A, 'fro') * norm(V, 'fro')));
