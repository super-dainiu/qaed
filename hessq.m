function [Q, H] = hessq(A)
% HESSQ Hessenberg form of a quaternion matrix.
%
% [Q, H] = HESSQ(A) computes the Hessenberg form of the quaternion matrix 
% A = A_0 + A_1*i + A_2*j + A_3*k where Q is a unitary quaternion matrix
% and H is an upper Hessenberg quaternion matrix such that: 
%
%       Q' * A * Q = H
%
% [Q, H] = HESSQ(A, RTOL) uses tolerance RTOL instead of the default eps.
%
% References:
%
% [1] Golub, G. H. and Van Loan, C. F., Matrix Computations, Algorithm
%     7.4.2, Johns Hopkins University Press, 1996. 
% 

narginchk(1, 1), nargoutchk(0, 2)

% Helper functions
householder_lapply = @(u, x) x - u * (u' * x);
householder_rapply = @(u, x) x - (x * u) * u';

n = size(A, 1);
Q = eyeq(n);

if ~isHessenberg(A)
    for i = 1:n-2
        v = A(i+1:end, i);
        r = eye(size(v));
        [u, ~] = householder_vector(v, r);
        
        % Apply to A
        A(i+1:end, i  :end) = householder_lapply(u, A(i+1:end, i  :end));
        A(   :   , i+1:end) = householder_rapply(u, A(   :   , i+1:end));
    
        % Apply to Q
        Q(   :   , i+1:end) = householder_rapply(u, Q(   :   , i+1:end));
    end
end

% Extract Hessenberg matrix
H = triu(A, -1);

if nargout < 2
    Q = H;
end

end

% Helper functions
function y = isHessenberg(A)
  % Returns true if A is upper Hessenberg
  y = all(abs(triu(A, -2)) < eps); 
end