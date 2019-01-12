function value = diffSpline(i, order, cfg, t)
d = cfg.k - 1;
len = cfg.nc;
if i > cfg.nc || i < 1
    error('Error. i should between [%d,%d], but input is %d\n', 1, cfg.nc, i);
end
r1 = cfg.knot(i+d) - cfg.knot(i);
r2 = cfg.knot(i+1+d) - cfg.knot(i+1);
if order == 1
    s1 = SplineBase(i, cfg.k-1, t, cfg.knot);
    if r1 == 0
        r1 = inf;
    end
    if r2 == 0
        r2 = inf;
    end
    if i == len
        value = d * (1/r1) * s1;
    else 
        s2 = SplineBase(i+1, cfg.k - 1, t, cfg.knot);
        value = d * ((1/r1) * s1 - (1/r2) * s2);
    end
else
    new_cfg = cfg;
    new_cfg.k = cfg.k - 1;
    s1 = diffSpline(i, order - 1, new_cfg, t);
    if r1 == 0
        r1 = inf;
    end
    if r2 == 0
        r2 = inf;
    end
    if i == len
        value = d * (1/r1) * s1;
    else 
        s2 = diffSpline(i+1, order - 1, new_cfg, t);
        value = d * ((1/r1) * s1 - (1/r2) * s2);
    end
end
