function [f, g] = cross(x1, y1, x2, y2)
f = x1*y2 - x2*y1;
g = zeros(4, 1);
g(1) = y2;
g(2) = -y1;
g(3) = -x2;
g(4) = x1;