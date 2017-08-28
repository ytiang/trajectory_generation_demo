tau = 0:0.001:1;
B = getSplineMatrix(1, tau, cfg);
px = p(1:cfg.nc);
py = p(cfg.nc+1:2*cfg.nc);
% B-spline curve
X = B * px;
Y = B * py;

% control polygon
Np = length(cfg.bpt) - 1;
for i = 1 : Np
    id0 = (i-1)*cfg.m;
    x = p(id0+1 : id0+cfg.k);
    y = p(id0+cfg.nc+1:id0+cfg.nc+cfg.k);
    if checkCoLinear(x, y) == 1
        index = [1:cfg.k, 1]';
    else
        index = convhull(x, y);
    end
    plot(x(index), y(index), 'g');
    hold on;
end
plot(px, py, 'm*', X, Y, 'k');
hold on;
for i = 1:4
    viscircles(region(i, 1:2), region(i, 3));
end
axis equal;
hold on;
% new_x = zeros(length(tau), 1);
% new_y = zeros(length(tau), 1);
% wx = w(1 : cfg.nc);
% wy = w(cfg.nc+1 : 2*cfg.nc);
% for i = 1 : length(tau)
%     paramx = nonUniformSplineParam(0, B(i, :), wx);
%     paramy = nonUniformSplineParam(0, B(i, :), wy);
%     new_x(i) = paramx.R * px;
%     new_y(i) = paramy.R * py;
% end
% plot(new_x, new_y, 'y');

