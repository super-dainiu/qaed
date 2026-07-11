function A = load_qmat(path)
% LOAD_QMAT Load a .qmat binary file (written by cpp/qaed_bench --save-out or
% scripts/save_qmat.m) as a qtfm quaternion matrix.
f = fopen(path, 'rb');
assert(f > 0, 'cannot open %s', path);
m = fread(f, 1, 'int64');
n = fread(f, 1, 'int64');
w = reshape(fread(f, m * n, 'double'), m, n);
x = reshape(fread(f, m * n, 'double'), m, n);
y = reshape(fread(f, m * n, 'double'), m, n);
z = reshape(fread(f, m * n, 'double'), m, n);
fclose(f);
A = quaternion(w, x, y, z);
end
