function group = groupByValue(A)
% GROUPBYVALUE Cluster indices for ordschur: equal values (within eps) share
% a group; groups numbered in ascending order of value.
[~, idx] = sort(A);
group = zeros(size(A));
current = 1;
for i = 2:length(A)
    if abs(A(idx(i)) - A(idx(i-1))) < eps
        group(idx(i)) = current;
    else
        current = current + 1;
        group(idx(i)) = current;
    end
end
group(idx(1)) = 1;
end
