function diff = getDSplineMatrix(order, T, cfg)
if order == 1
    diff = zeros(length(T), (cfg.nc-1));
    for i = 1:length(T)
        for j = 1:cfg.nc-1
            coe = SplineBase(j, cfg.k-1, T(i), cfg.knot(2:length(cfg.knot)-1));
            diff(i, j) = (cfg.k-1)*coe/(cfg.knot(j+cfg.k) - cfg.knot(j+1));
        end
    end
end