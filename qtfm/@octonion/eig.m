function [V, D] = eig(A) %#ok<STOUT,INUSD>
% EIG    Eigenvalues and eigenvectors.
% (Octonion overloading of standard Matlab function.)
%
% Acceptable calling sequences are: [V,D] = EIG(X) and V = EIG(X).
% The results are as for the standard Matlab EIG function (q.v.).

% Copyright © 2020 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

unimplemented(mfilename)

% TODO See the following paper, which although it studies eigenvalues and
% eigenvectors of octonion matrices (small ones!), it doesn't give an
% algorithm for computing them:
%
% Dray, Tevian and Manogue, Corinne, A. (1998),
% 'Finding octonionic eigenvalues using Mathematica',
% Computer Physics Communications, 115, 536–547.
% DOI :10.1016/S0010-4655(98)00126-X
% Also available as arXiv:math/9807133v1 23 July 1998.

end

% $Id: eig.m 1113 2021-02-01 18:41:09Z sangwine $
