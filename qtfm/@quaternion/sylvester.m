function x = sylvester(a, b, c)
% SYLVESTER  Solve quaternion sylvester equation ax + xb = c for x.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2016, 2020 Stephen J. Sangwine.
% See the file : Copyright.m for further details.

% Reference:
%A
% Drahoslava Janovska and Gerhard Opfer,
% "Linear equations in quaternionic variables",
% Mitteilungen der Mathematischen Gesellschaft in Hamburg 27 (2008),
% 223–234. Theorem 2.3. Available: https://www.math.uni-hamburg.de/mathges/

narginchk(3, 3), nargoutchk(0, 1)

if ~isa(a, 'quaternion') || ~isa(b, 'quaternion') || ~isa(c, 'quaternion')
    error('All parameters must be quaternions.')
end

na = normq(a);
nb = normq(b);

if na <= nb
    % Use the second formula in Theorem 2.3 of the paper above.
    
    ib = b.^-1;
    
    x = (c + conj(a) .* c .* ib) .* (2 .* s(a) + b + na .* ib).^-1;
    
else
    % Use the first formula.
    
    ia = a.^-1;
    
    x = (2 .* s(b) + a + nb .* ia).^-1 .* (c + ia .* c .* conj(b));
end

end

% Note: An earlier version of this function released in QTFM 2.5 used
% adjoint matrices and the MATLAB sylvester function.

% TODO This code works only for scalar parameters a, b, c, whereas the
% previous code based on adjoints would work for matrices of quaternions!
% Restore the previous code and limit the Janovska and Opfer code to the
% scalar case (for better accuracy and faster computation), but use the
% adjoint method and the Matlab sylvester function in order to support the
% matrix case. NB What if you want to solve elementwise?

% TODO See also: D. Shirokov, 'Basis-free Solution to Sylvester Equation in
% Clifford Algebra of Arbitrary Dimension', Advances in Applied Clifford
% Algebras 31, 70 (2021). DOI:10.1007/s00006-021-01173-0 [Section 1 covers
% the quaternion case.]

% $Id: sylvester.m 1131 2021-12-21 21:52:48Z sangwine $
