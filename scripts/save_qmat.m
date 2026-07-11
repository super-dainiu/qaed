function save_qmat(path, A)
% SAVE_QMAT Save a qtfm quaternion matrix in the .qmat binary format used by
% the C++ port (cpp/qaed_bench --load).
%
% Format: int64 m, int64 n, then column-major double arrays w, x, y, z.
A = quaternion(A);
f = fopen(path, 'wb');
assert(f > 0, 'cannot open %s', path);
[m, n] = size(A);
fwrite(f, int64(m), 'int64');
fwrite(f, int64(n), 'int64');
w = A.w; if isempty(w), w = zeros(m, n); end
fwrite(f, w, 'double');
fwrite(f, A.x, 'double');
fwrite(f, A.y, 'double');
fwrite(f, A.z, 'double');
fclose(f);
end
