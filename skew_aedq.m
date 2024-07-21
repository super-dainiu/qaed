function [Q, D] = skew_aedq(D, E, rtol, verbose, bar)
narginchk(2, 4), nargoutchk(0, 2)

if nargin < 5
    bar = true;
end%if

if nargin < 4
    verbose = false;
end%if

if nargin < 3
    rtol = eps;
end%if

H      = diag(D) - diag(conj(E), 1) + diag(E, -1);  A = H;
[n, ~] = size(D); ilo = 1; ihi = n; Q = eyeq(n);

Q_time = 0;

if bar
    b = waitbar(0, strcat('tridiag aed, n = ', num2str(n)));
end%if

while ihi > 1
    while ihi - ilo > 2
        
        if bar
            waitbar(1 - ihi / n, b)
        end%if

        for j = ihi - 1 : -1 : ilo + 1
            if abs(H(j, j - 1)) <= rtol * (abs(H(j, j)) + abs(H(j - 1, j - 1)))
                ilo = j;
                break
            end%if
        end%for

        ns  = aed_num_shifts(ihi - ilo);
        win  = aed_win_size(ihi - ilo, ns);
        whi = ihi; wlo = max(ihi - win + 1, ilo); sp  = wlo - 1;

        if ihi - ilo + 1 >= aed_min_size() && sp > ilo && win > 4
            [U, D] = tridiagq(diag(H(wlo : whi, wlo : whi)), diag(H(wlo : whi, wlo : whi), -1), rtol, false, false);
            spike  = U' * H(wlo : whi, sp); 
            [~, index] = sort(abs(spike), "descend"); spike = spike(index);
            H(wlo : whi, wlo : whi) = diag(D(index));

            T0 = tic;
            Q(  :  , wlo : whi) = Q(  :  , wlo : whi) * U;
            Q(  :  , wlo : whi) = Q(  :  , sp + index);
            Q_time = Q_time + toc(T0);

            for i = win : -1 : 1
                if abs(spike(i)) <= rtol * abs(D(i))
                    spike(i) = 0; whi = whi - 1;
                end%if
            end%for
            H(wlo : ihi, sp) = spike; H(sp, wlo : ihi) = -spike';
            shifts = diag(H(wlo : whi, wlo : whi)); ihi = whi;
            shifts = shifts(end : -1 : 1);

            % hessenberg reduction
            [U, H(sp : ihi, sp : ihi)] = hessq(H(sp : ihi, sp : ihi), rtol);
            T0 = tic;
            Q(  :  , sp : ihi) = Q(  :  , sp : ihi) * U;
            Q_time = Q_time + toc(T0);

            if numel(shifts) / win < (1 - nibble())   % if sufficient deflation
                continue
            end%if
        else
            shifts = [H(ihi, ihi)];
        end%if

        %% case 3: implicit qr steps

        ls = 0;
        ns = min(numel(shifts), ns);

        while ihi - ilo > 2 && ls < ns
            % small-bulge shifts
            ls     = ls + 1;
            x      = shifts(ls);
            si     = - 2 * x.w;
            ti     = abs(x) ^ 2;
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
            
            while ihi - ilo > 2 && ...
                    (abs(H(ihi, ihi - 1)) <= ...
                    rtol * (abs(H(ihi, ihi)) + abs(H(ihi - 1, ihi - 1))))
                H(ihi, ihi - 1) = 0;
                ihi  =  ihi - 1;
            end%while
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

%--------------------------------------------------------------------------

function NS = aed_num_shifts(ihi)
if ihi < 30
    NS = 2;
elseif ihi < 60
    NS = 4;
elseif ihi < 150
    NS = 10;
elseif ihi < 590
    NS = max(0, ihi / round(log2(double(ihi))));
elseif ihi < 3000
    NS = 64;
else
    NS = 128;
end%if

NS = max(2, NS - mod(NS, 2));
end%function

%--------------------------------------------------------------------------

function WS = aed_win_size(ihi, NS)
if ihi <= 500
    WS = NS;
else
    WS = (3 * NS) / 2;
end%if

WS = max(4, WS - mod(WS, 2));
end%function

%--------------------------------------------------------------------------

function MS = aed_min_size()
MS = 12;
end%function

%--------------------------------------------------------------------------

function N = nibble()
N = 0.14;
end%function


    

