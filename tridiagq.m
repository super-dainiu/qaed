function [Q, D] = tridiagq(D, E, rtol, verbose, bar)
narginchk(2, 5), nargoutchk(0, 2)

if nargin < 5
    bar = true;
end%if

if nargin < 4
    verbose = false;
end%if

if nargin < 3
    rtol = eps;
end%if

H      = diag(D) - diag(conj(E), 1) + diag(E, -1);
[n, ~] = size(D); ilo = 1; ihi = n; Q = eyeq(n);
Q_time = 0;

if bar
    b = waitbar(0, strcat('tridiag, n = ', num2str(n)));
end%if

while ihi > 1
    while ihi - ilo > 2

        for j = ihi - 1 : -1 : ilo + 1
            if abs(H(j, j - 1)) <= rtol * (abs(H(j, j)) + abs(H(j - 1, j - 1)))
                ilo = j;
                break
            end%if
        end%for
    
        x  = H(ihi, ihi);
        si = - 2 * x.w;
        ti = abs(x) ^ 2;
        v      = [
            H(ilo, ilo) ^ 2  + H(ilo, ilo + 1) * H(ilo + 1, ilo) + si * H(ilo, ilo) + ti; ...
            H(ilo + 1, ilo) * H(ilo, ilo) + H(ilo + 1, ilo + 1) * H(ilo + 1, ilo) + si * H(ilo + 1, ilo); ...
            H(ilo + 2, ilo + 1) * H(ilo + 1, ilo);
        ];
        r      = zeros(numel(v), 1); r(1) = 1;
        [u, ~] = householder_vector(v, r);
    
        H(ilo : ilo + 3, ilo : ilo + 2) = H(ilo : ilo + 3, ilo : ilo + 2) - (H(ilo : ilo + 3, ilo : ilo + 2) * u) * u';
        H(ilo : ilo + 2, ilo : ilo + 3) = H(ilo : ilo + 2, ilo : ilo + 3) - u * (u' * H(ilo : ilo + 2, ilo : ilo + 3));
        
        T0 = tic;
        Q(  :  , ilo : ilo + 2) = Q(  :  , ilo : ilo + 2) - (Q(  :  , ilo : ilo + 2) * u) * u';
        Q_time = Q_time + toc(T0);

        for i = ilo : (ihi - 2)
            e      = min(i + 3, ihi);
            ee     = min(e + 1, ihi);
            v      = H(i + 1 : e, i);
            r      = zeros(numel(v), 1); r(1) = 1;
            [u, ~] = householder_vector(v, r);
    
            H( i + 1 : e , i : ee) = H( i + 1 : e , i : ee) - u * ...
                                    (u' * H( i + 1 : e , i : ee));
            H( i : ee , i + 1 : e) = H( i : ee , i + 1 : e) - ...
                                    (H( i : ee , i + 1 : e) * u) * u';
            T0 = tic;
            Q(   :    , i + 1 : e) = Q(   :    , i + 1 : e) - ...
                                    (Q(   :    , i + 1 : e) * u) * u';
            Q_time = Q_time + toc(T0);
        end%for
        
        while ihi > ilo + 1 && ...
                (abs(H(ihi, ihi - 1)) <= ...
                rtol * (abs(H(ihi, ihi)) + abs(H(ihi - 1, ihi - 1))))
            H(ihi, ihi - 1) = 0;
            ihi  =  ihi - 1;
            if bar
                waitbar(1 - ihi / n, b)
            end%if
        end%while
    
        H = triu(H, -1);
        H = tril(H,  1);
    end%while
    [U, H( ilo : ihi , ilo : ihi )] = iqrq(H( ilo : ihi , ilo : ihi ), rtol, false, false, false);
    T0 = tic;
    Q(  :  , ilo : ihi ) = Q( : , ilo : ihi ) * U;
    Q_time = Q_time + toc(T0);
    ihi = ilo - 1; ilo = 1;
end%while

if bar
    delete(b);
end%if

if verbose
    fprintf('time to construct Q = %f s\n', Q_time);
end%if

if nargout < 2
    Q = diag(H);
else
    D = diag(H);
end%if
end%function
    

