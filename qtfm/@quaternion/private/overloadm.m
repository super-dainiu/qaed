function R = overloadm(F, Q)
% Private function to implement overloading of matrix Matlab functions.
% Called to apply the function F to the quaternion array Q by operating on
% an adjoint representation of Q. F must be a string, giving the name of
% the function F. The calling function can pass this string using
% mfilename, for simplicity of coding.

% Copyright © 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

H = str2func(F); % A handle to the function designated by F.

R = unadjoint(H(adjoint(Q, 'real')), 'real');

end

% Implementation note: the choice of the real adjoint is to avoid returning
% a complex result in the case where Q is a real quaternion. Because the
% adjoint method can be inaccurate, the result of the unadjoint may have
% small non-zero imaginary parts, resulting in a biquaternion result from a
% quaternion matrix.

% Implementation note: the use of the 'block' option to un/adjoint was
% tried Jan 2011, but results showed that although it sometimes gave a
% slightly better result, on average it didn't.

% $Id: overloadm.m 1113 2021-02-01 18:41:09Z sangwine $
