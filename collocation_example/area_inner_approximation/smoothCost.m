% function [f, g] = smoothCost(w, p, B, cfg)
function [f] = smoothCost(w, p, B, cfg)
wx = w(1 : cfg.nc);
wy = w(cfg.nc+1 : 2*cfg.nc);
px = p(1 : cfg.nc);
py = p(cfg.nc+1 : 2*cfg.nc);
[rows, cols] = size(B);
f = 0;
g = zeros(cfg.nc*2, 1);
for i = 1 : rows - 2
    paramx1 = nonUniformSplineParam(0, B(i, :), wx);
    paramx2 = nonUniformSplineParam(0, B(i+1, :), wx);
    paramx3 = nonUniformSplineParam(0, B(i+2, :), wx);
    paramy1 = nonUniformSplineParam(0, B(i, :), wy);
    paramy2 = nonUniformSplineParam(0, B(i+1, :), wy);
    paramy3 = nonUniformSplineParam(0, B(i+2, :), wy);
    x1 = paramx1.R * px;
    x2 = paramx2.R * px;
    x3 = paramx3.R * px;
    y1 = paramy1.R * py;
    y2 = paramy2.R * py;
    y3 = paramy3.R * py;
    dp1 = [x2-x1, y2-y1];
    dp2 = [x3-x2, y3-y2];
    ddp = dp2 - dp1;
    f = f + ddp * ddp';
    g(1:cfg.nc) = g(1:cfg.nc) + ...
        2 * ddp(1) * (paramx3.W - 2*paramx2.W + paramx1.W) * px;
    g(cfg.nc+1 : 2*cfg.nc) = g(cfg.nc+1 : 2*cfg.nc) + ...
        2 * ddp(2) * (paramy3.W - 2*paramy2.W + paramy1.W) * py;
end