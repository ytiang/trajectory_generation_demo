function [c, ceq, dc, dceq] = nlonConstraints(params, cfg, region)
Np = length(cfg.bpt) - 1;
p_num = length(params) / 2;
multi = cfg.k - cfg.m;
c = zeros(Np*(cfg.k + 1) + (Np-2)*2*(cfg.k-2*multi)+2*cfg.m, 1);
dc = zeros(Np*(cfg.k + 1) + (Np-2)*2*(cfg.k-2*multi)+2*cfg.m, 2*p_num);
ceq = [];
dceq = [];
for i = 1 : Np
    id0 = (i-1)*cfg.m;
    x = params(id0+1 : id0+cfg.k);
    y = params(id0+p_num+1:id0+p_num+cfg.k);
    circle = region(i, :);
    for j = 1 : cfg.k
        c((i-1)*cfg.k + j) = (x(j) - circle(1))^2 + (y(j) - circle(2))^2 - circle(3)^2;
        dc((i-1)*cfg.k + j, id0 + j) = 2 * (x(j)-circle(1));
        dc((i-1)*cfg.k + j, id0 + p_num + j) = 2 * (y(j)-circle(2));
    end
    if checkCoLinear(x, y) == 1
        index = [1:cfg.k, 1]';
    else
        index = convhull(x, y);
    end
    [area, dx, dy] = polyArea(x, y, index);
    c(Np*cfg.k + i) = -area;
    dc(Np*cfg.k + i, id0+1 : id0+cfg.k) = -dx';
    dc(Np*cfg.k + i, id0+p_num+1 : id0+p_num+cfg.k) = -dy';
    if i == 1
        for j = 1 : cfg.m
            c(Np*(cfg.k + 1)+j) = ... 
                region(i+1, 3)^2 - (x(j)-region(i+1, 1))^2 - (y(j) - region(i+1, 2))^2;
            dc(Np*(cfg.k + 1)+j, id0+j) = -2*(x(j)-region(i+1, 1));
            dc(Np*(cfg.k + 1)+j, id0+p_num+j) = -2*(y(j)-region(i+1, 2));
        end
    elseif i == Np
        for j = multi+1 : cfg.k
            c(Np*(cfg.k + 1) + (Np-2)*2*(cfg.k-2*multi)+cfg.m+j) = ... 
                region(i-1, 3)^2 - (x(j)-region(i-1, 1))^2 - (y(j) - region(i-1, 2))^2;
            dc(Np*(cfg.k + 1) + (Np-2)*2*(cfg.k-2*multi)+cfg.m+j, id0+j) = ... 
                -2*(x(j)-region(i-1, 1));
            dc(Np*(cfg.k + 1) + (Np-2)*2*(cfg.k-2*multi)+cfg.m+j, id0+p_num+j) = ... 
                -2*(y(j)-region(i-1, 2));
        end
    else
        for j = multi+1 : cfg.k - multi
            c(Np*(cfg.k + 1) + cfg.m + j) = ... 
                region(i-1, 3)^2 - (x(j)-region(i-1, 1))^2 - (y(j) - region(i-1, 2))^2;
            dc(Np*(cfg.k + 1) + cfg.m + j, id0+j) = ... 
                -2*(x(j)-region(i-1, 1));
            dc(Np*(cfg.k + 1) + cfg.m + j, id0+p_num+j) = ... 
                -2*(y(j)-region(i-1, 2));
            c(Np*(cfg.k + 1) + (Np-2)*(cfg.k-2*multi)+cfg.m+j) = ... 
                region(i+1, 3)^2 - (x(j)-region(i+1, 1))^2 - (y(j) - region(i+1, 2))^2;
            dc(Np*(cfg.k + 1) + (Np-2)*(cfg.k-2*multi)+cfg.m+j, id0+j) = ... 
                -2*(x(j)-region(i+1, 1));
            dc(Np*(cfg.k + 1) + (Np-2)*(cfg.k-2*multi)+cfg.m+j, id0+p_num+j) = ... 
                -2*(y(j)-region(i+1, 2));
        end
    end
end
