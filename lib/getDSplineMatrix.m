function M = getDSplineMatrix(order, dim, T, cfg)
if order == 1
    M = zeros(length(T)*dim, (cfg.nc-1)*dim);
    for i = 1:length(T)
        for j = 1:cfg.nc-1
            if(i == length(T) && j == cfg.nc - 1)
                coe = 1;
            else
                coe = SplineBase(j, cfg.k-1, T(i), cfg.knot(2:length(cfg.knot)-1));
            end
            coe = (cfg.k-1)*coe/(cfg.knot(j+cfg.k) - cfg.knot(j+1));
            M((i-1)*dim+1:i*dim, (j-1)*dim+1:j*dim) =coe*eye(dim); 
        end
    end
end