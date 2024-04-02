function x = sylvesterc(a, b, c)
% SYLVESTERC Solver for a Special Case of the Sylvester Equation for
% Quaternion Scalars
%
% x = SYLVESTERC(a, b, c) solves a specialized Sylvester equation
%
%   ax - xb = c
%
% where a and b are complex scalars, and c is a quaternion scalar. The
% function computes the solution x, which is a quaternion scalar. 
%
% References:
%
% [1] Zhang, Fuzhen. "Quaternions and matrices of quaternions.” 
%     Linear Algebra and its Applications 251 (1997): 21-57.

narginchk(3, 3), nargoutchk(1, 1)

a_1 = a.s + 1i * a.x;
b_1 = b.s + 1i * b.x;
b_2 = b.s - 1i * b.x;
c_1 = c.s + 1i * c.x;
c_2 = c.y + 1i * c.z;

x_1 = c_1 / (a_1 - b_1);
x_2 = c_2 / (a_1 - b_2);
x   = quaternion(real(x_1), imag(x_1), real(x_2), imag(x_2));
end