function M = getSplineMatrix(dim, T, cfg)
row = length(T);
M = zeros(row*dim, cfg.nc*dim);
for i=1:row
    for j=1:cfg.nc
        if(j == cfg.nc && i == row)
            coe = 1;
        else
            coe = SplineBase(j, cfg.k, T(i), cfg.knot);
        end
        M((i-1)*dim+1:i*dim, (j-1)*dim+1:j*dim) = coe * eye(dim);
    end
end