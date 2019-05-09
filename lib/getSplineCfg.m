function cfg = getSplineCfg(k, seg_num, multi)
cfg.k = k;
cfg.m = multi;
bpts = zeros(seg_num+1, 1);
knot = [];
for i=1:seg_num+1
    bpts(i) = 1 / (seg_num) * (i-1);
    if(i==1 || i==seg_num+1)
        m = k;
    else
        m = multi;
    end
    knot = [knot; bpts(i)*ones(m, 1)];
end
cfg.knot = knot;
cfg.bpt = bpts;
cfg.nc = length(knot) - k;
