function [w, x, y, z] = qmat_parts(A)
% QMAT_PARTS Extract the four real component arrays of a qtfm quaternion
% matrix (helper for the *_cpp wrappers). Pure quaternions get w = zeros.
if ~isa(A, 'quaternion')
    A = quaternion(A);
end
x = A.x; y = A.y; z = A.z;
w = A.w;
if isempty(w), w = zeros(size(x)); end
end
