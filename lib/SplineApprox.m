function value = SplineApprox(ctrp, dim, t, cfg)
n = length(ctrp) / dim;
M = zeros(dim, n*dim);
value = zeros(length(t), dim);
pt = zeros(dim, 1);
for j=1:length(t)
    for i=1:n
        if(t(j) == cfg.knot(length(cfg.knot)) && i == n)
            coe = 1;
        else
            coe = SplineBase(i, cfg.k, t(j), cfg.knot);
        end
        M(:, (i-1)*dim+1:i*dim) = coe*eye(dim);
    end
    pt = M * ctrp;
    value(j,:) = pt';
end