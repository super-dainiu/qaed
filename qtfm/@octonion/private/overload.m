function R = overload(F, Q, varargin)
% Private function to implement overloading of Matlab functions. Called to
% apply the function F to the octonion array Q by operating on the
% quaternion components of Q with F. F must be a string, giving the name of
% the function F. The calling function can pass this string using
% mfilename, for simplicity of coding. The varargin parameters must not be
% octonions, since this code doesn't know whether to split them into
% quaternions or pass them unchanged.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

H = str2func(F); % A handle to the function designated by F.

R = Q;

R.a = H(R.a, varargin{:});
R.b = H(R.b, varargin{:});

end

% $Id: overload.m 1113 2021-02-01 18:41:09Z sangwine $
