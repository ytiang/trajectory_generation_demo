function M = differentialMatrix(order, n, dim)
if(order == 1)
    M =zeros((n-1)*dim,n*dim);
    for i=1:n-1
        M((i-1)*dim+1:i*dim, (i-1)*dim+1:(i+1)*dim) = [-eye(dim), eye(dim)];
    end
end