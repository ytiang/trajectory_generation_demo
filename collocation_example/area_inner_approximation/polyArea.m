function [area, dx, dy] = polyArea(x, y, index)
num = length(x);
area = 0.0;
dx = zeros(num, 1);
dy = zeros(num, 1);
for i = 1 : length(index)-1
    [v, dp] = cross(x(index(i)), y(index(i)), x(index(i+1)), y(index(i+1)));
    area = area + v;
    dx(index(i)) = dx(index(i)) + dp(1);
    dx(index(i+1)) = dx(index(i+1)) + dp(2);
    dy(index(i)) = dy(index(i)) + dp(3);
    dy(index(i+1)) = dy(index(i+1)) + dp(4);
end
area = area / 2;
dx = dx / 2;
dy = dy / 2;