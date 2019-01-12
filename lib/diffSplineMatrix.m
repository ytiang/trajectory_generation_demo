function dspline = diffSplineMatrix(order, cfg, t)
dspline = zeros(length(t), cfg.nc);
for i = 1:length(t)
    for j = 1 : cfg.nc
        dspline(i, j) = diffSpline(j, order, cfg, t(i));
    end    
end