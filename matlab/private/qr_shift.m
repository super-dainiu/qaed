function s = qr_shift(A)
% QR_SHIFT Wilkinson-type shift from the trailing 2x2 block: the standardized
% eigenvalue class of A closest to the standardized (2,2) entry.
[~, T1] = schur1(A(2, 2));
[~, T2] = schur2(A);
if abs(T2(1, 1) - T1) < abs(T2(2, 2) - T1)
    s = T2(1, 1);
else
    s = T2(2, 2);
end
end
