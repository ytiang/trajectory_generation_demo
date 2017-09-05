function block = Query(A, i, j, n)
if j == 1
    block = A(i, 1:n(1));
else
    block = A(i, sum(n(1:j-1))+1 : sum(n(1:j)));
end