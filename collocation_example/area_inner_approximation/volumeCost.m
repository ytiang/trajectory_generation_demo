function [f, g] = volumeCost(params, cfg)
Np = length(cfg.bpt) - 1;
p_num = length(params) / 2;
f = 0.0;
g = zeros(p_num*2, 1);
for i = 1 : Np
    id0 = (i-1)*cfg.m;
    x = params(id0+1 : id0+cfg.k);
    y = params(id0+p_num+1:id0+p_num+cfg.k);
    if checkCoLinear(x, y) == 1
        index = [1:cfg.k, 1]';
    else
        index = convhull(x, y);
    end
    % the 2-volume of convex hull
    [area, dx, dy] = polyArea(x, y, index);
    f = f + area;
    g(id0+1 : id0+cfg.k) = g(id0+1 : id0+cfg.k) + dx;
    g(id0+p_num+1:id0+p_num+cfg.k) = g(id0+p_num+1:id0+p_num+cfg.k) + dy;
end
f = -f;
g = -g;