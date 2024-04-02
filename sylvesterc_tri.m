function X = sylvesterc_tri(T, b, c)
% SYLVESTERC_TRI Solver for Sylvester Equation with Triangular
% Coefficient
%
% X = SYLVESTERC_TRI(T, b, c) solves the Sylvester equation of the form:
% TX - Xb = c, where T is an upper triangular matrix with complex diagonal
% elements. The function is only supported for cases where T is upper
% triangular.
%
% Example Usage:
% T = triu(randq(4)); % Example upper triangular matrix
% b = randq(1); % Example scalar
% c = randq(4, 1); % Example vector
% x = sylvesterc_tri(T, b, c); % Solving the Sylvester equation
%
% References:
%
% [1] Zhang, Fuzhen. "Quaternions and matrices of quaternions.” 
%     Linear Algebra and its Applications 251 (1997): 21-57.

narginchk(3, 3), nargoutchk(1, 1)

n = size(T, 1); X = zerosq(n, 1);

for i = n:-1:1
    X(i) = sylvesterc(T(i, i), b, c(i));
    c = c(1:i-1) - T(1:i-1, i) * X(i);
end
end

