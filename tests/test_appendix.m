% TEST_APPENDIX Timing comparison of the three implementations for the
% appendix tables: C++ via MEX (dispatch on) and pure MATLAB (dispatch off),
% on the same matrices as the main-text C++ experiments.
%
% Prints one machine-parsable line per run:
%   APPENDIX type=<full|hess> n=<n> alg=<aed|iqr> impl=<mex|matlab> time=<s>

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab'), fullfile(root, 'qtfm'), fullfile(root, 'cpp'));

assert(exist('qaed_mex', 'file') == 3, 'qaed_mex not built / not on path');
rtol = eps;

for type = ["full", "hess"]
    % MEX-backed MATLAB: same range as the main-text tables (up to 4096).
    for n = [64, 128, 256, 512, 1024, 2048, 4096]
        H = make_matrix(type, n);
        qaed_accel(true);
        tic; [~, ~] = aedq(H, rtol); t = toc;
        fprintf('APPENDIX type=%s n=%d alg=aed impl=mex time=%.3f\n', type, n, t);
        if n <= 2048
            tic; [~, ~] = iqrq(H, rtol); t = toc;
            fprintf('APPENDIX type=%s n=%d alg=iqr impl=mex time=%.3f\n', type, n, t);
        end
    end
    % Pure MATLAB: up to 1024 (larger sizes are impractical).
    for n = [64, 128, 256, 512, 1024]
        H = make_matrix(type, n);
        qaed_accel(false);
        tic; [~, ~] = aedq(H, rtol); t = toc;
        fprintf('APPENDIX type=%s n=%d alg=aed impl=matlab time=%.3f\n', type, n, t);
        tic; [~, ~] = iqrq(H, rtol); t = toc;
        fprintf('APPENDIX type=%s n=%d alg=iqr impl=matlab time=%.3f\n', type, n, t);
    end
end
qaed_accel(true);

function H = make_matrix(type, n)
rng(0);
A = randq(n) .* randn(n);
if type == "hess"
    H = triu(A, -1);
else
    qaed_accel(true);
    [~, H] = hessq(A);
end
end
