function test_qfft
% Test code for the fast quaternion Fourier transform.

% This code tests the following functions:
%
%  qfft     qfft2
% iqfft    iqfft2
%  qfft2ds   fft2ds
% iqfft2ds  ifft2ds
% iqdft2ds  qdft2ds

% Copyright © 2005, 2010, 2022 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% We have to test all of the above functions separately, since they do not
% call each other, as is the case with the corresponding dqft functions. In
% addition, we need to verify the results of the FFTs against the DFT code
% in the corresponding qdft functions, e.g. for qfft2, this is qdft2. For
% each transform we need to verify that it inverts correctly for all four
% combinations of real/complex data and real/complex axis, and for left and
% right exponentials, and we need to verify each against the corresponding
% DFT code. We do this with the default axes defined in the private
% function dft_axis, and with random axes, to check that the decompositions
% used in the qfft code are correct.

disp('Testing single-sided fast quaternion Fourier transforms ...')

% Define one real and one complex quaternion array.

q = randq(10) .* randn(10);
b = complex(randq(10), randq(10)) .* randn(10);

T = 1e-10;

% Construct a set of test axes, some fixed, some random.

RA = [unit(quaternion(1,1,1)), qi, qj, qk, randv, randv, randv, randv];
CA = [complex(quaternion(1,1,1), quaternion(0,1,-1)), ...
      unit(complex(quaternion(0,1,1), quaternion(1,1,1))), ...
      unit(qi + qj + qk .* 1i), unit(qi + randv .* 1i), ...
      unit(complex(randv, randv)), unit(complex(randv, randv)), ...
      unit(complex(randv, randv)), unit(complex(randv, randv))];
  
assert(length(RA) == length(CA));

if  isreal(CA); error('Complex axis is not complex.'); end
if ~isreal(RA); error('Real axis is complex.'); end

for k = 1:length(RA)
    
    R = RA(k);
    C = CA(k);
    
    % 1D FFT code against its own inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqfft(qfft(q, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 1L with axis: ', char(R)]);
    compare(q, iqfft(qfft(q, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 1R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqfft(qfft(q, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 2L with axis: ', char(R)]);
    compare(q, iqfft(qfft(q, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 2R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqfft(qfft(b, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 3L with axis: ', char(R)]);
    compare(b, iqfft(qfft(b, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 3R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqfft(qfft(b, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 4L with axis: ', char(R)]);
    compare(b, iqfft(qfft(b, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 4R with axis: ', char(R)]);
    
    % 2D FFT code against its own inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqfft2(qfft2(q, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 1L with axis: ', char(R)]);
    compare(q, iqfft2(qfft2(q, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 1R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqfft2(qfft2(q, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 2L with axis: ', char(R)]);
    compare(q, iqfft2(qfft2(q, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 2R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqfft2(qfft2(b, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 3L with axis: ', char(R)]);
    compare(b, iqfft2(qfft2(b, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 3R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqfft2(qfft2(b, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 4L with axis: ', char(R)]);
    compare(b, iqfft2(qfft2(b, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 4R with axis: ', char(R)]);
    
    % 1D FFT code against DFT inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqdft(qfft(q, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 5L with axis: ', char(R)]);
    compare(q, iqdft(qfft(q, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 5R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqdft(qfft(q, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 6L with axis: ', char(R)]);
    compare(q, iqdft(qfft(q, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 6R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqdft(qfft(b, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 7L with axis: ', char(R)]);
    compare(b, iqdft(qfft(b, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 7R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqdft(qfft(b, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 8L with axis: ', char(R)]);
    compare(b, iqdft(qfft(b, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 8R with axis: ', char(R)]);
    
    % 2D FFT code against DFT inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqdft2(qfft2(q, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 5L with axis: ', char(R)]);
    compare(q, iqdft2(qfft2(q, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 5R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqdft2(qfft2(q, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 6L with axis: ', char(R)]);
    compare(q, iqdft2(qfft2(q, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 6R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqdft2(qfft2(b, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 7L with axis: ', char(R)]);
    compare(b, iqdft2(qfft2(b, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 7R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqdft2(qfft2(b, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 8L with axis: ', char(R)]);
    compare(b, iqdft2(qfft2(b, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 8R with axis: ', char(R)]);

end

disp('Passed');

% Now test the double-sided transforms. We use a different approach here,
% in order to verify the results against a matrix algorithm published in:
%
% Sangwine, S. J. and Ell, T. A., ‘Complex and Hypercomplex Discrete
% Fourier Transforms Based on Matrix Exponential Form of Euler's Formula’,
% Applied Mathematics and Computation, 219(2), October 2012, 644–655.
% doi:10.1016/j.amc.2012.06.055.

disp('Testing double-sided fast quaternion Fourier transforms ...')

% Test with real data from the array q, first checking inversion, then
% checking against the matrix code.

A = {randv, randv}; % Transform axes.

compare(q, iqfft2ds(qfft2ds(q, A), A), T, 'qfft2ds failed test 9');
compare(q,  ifft2ds( fft2ds(q)),       T, 'fft2ds failed test 10');

p = qfft2ds(q, A);
Q = adjoint(q, 'real', 'block');
P = matdft2(Q, adjoint(A{1}, 'real'), adjoint(A{2}, 'real'));

compare(adjoint(p, 'real', 'block'), P, T, 'qfft2ds failed test 11');

% Test with complex data from the array b.

B = cellfun(@unit, {complex(randv, randv), complex(randv, randv)}, ...
    'UniformOutput', false); % Transform axes.

compare(b, iqfft2ds(qfft2ds(b, B), B), T, 'qfft2ds failed test 12');
compare(b,  ifft2ds( fft2ds(b)),       T,  'fft2ds failed test 13');

p = qfft2ds(b, B);
Q = adjoint(b, 'real', 'block');
P = matdft2(Q, adjoint(B{1}, 'real'), adjoint(B{2}, 'real'));

compare(adjoint(p, 'real', 'block'), P, T, 'qfft2ds failed test 14')

% Now check against the DFT code. We use a small piece of q for this, as
% the DFT code is slow.

q = q(1:3, 1:3);

compare(qdft2ds(q, A), qfft2ds(q, A), T, 'qdftd2s/qfft2ds failed test 15');
compare(q, iqdft2ds(qdft2ds(q, A), A), T, 'qdft2ds/iqdft2ds failed test 16');

disp('Passed');

end

% Quaternion two-sided DFT based on matrix exponential, used in the last
% few tests above to provide an assurance that the QDFT2DS code is
% correctly implementing the two-sided quaternion DFT.

function F = matdft2(f, J, K)
% This function is a matrix implementation of a two-sided two-dimensional
% Fourier transform as Theorem 4.1 in Todd Ell's thesis. The input array f
% must be stored as a block matrix with blocks representing hypercomplex
% elements, e.g. 4x4 blocks for the quaternion case.

assert(all(size(J) == size(K))); % Check that J and K are square and that
assert(size(J,1) == size(J,2));  % they are the same size.

A = size(J, 1);      % This is the size of the matrix block element,
                     % representing one (hyper)complex sample.

M = size(f, 1) ./ A; % The number of (hyper)complex elements horizontally.
N = size(f, 2) ./ A; % The number of (hyper)complex elements vertically.

assert(M .* A == size(f, 1)); % Check that the dimensions of f in both
assert(N .* A == size(f, 2)); % directions are integral multiples of the
                              % size of J and K.
F = zeros(size(f), class(f));

for u = 0:M-1 %#ok<BDSCI> 
    for v = 0:N-1 %#ok<BDSCI> 
        for m = 0:M-1 %#ok<BDSCI> 
            for n = 0:N-1 %#ok<BDSCI> 
                % The four loop variables are the mathematical indices of
                % the elements of f and F. Each element is a A-by-A matrix.
                % Therefore, the indices to be used to slice f and F need
                % to be adjusted to be multiples of A, as well as indexing
                % from 1 rather than 0.
                
                F(A * u + 1:A * u + A, A * v + 1:A * v + A) = ...
                F(A * u + 1:A * u + A, A * v + 1:A * v + A) + ...
                Expm(-J, 2 .* pi .* m .* u ./ M) * ...
                f(A * m + 1:A * m + A, A * n + 1:A * n + A) * ...
                Expm(-K, 2 .* pi .* n .* v ./ N);
            end
        end
    end
end

end

function E = Expm(R, theta)
% Returns the matrix expm(R .* theta) using trig functions rather than
% expm. R must be a square root of -eye(size(R)) (not verified).

E = eye(size(R)) .* cos(theta) + R .* sin(theta);
end

% $Id: test_qfft.m 1168 2022-11-20 19:07:25Z sangwine $
