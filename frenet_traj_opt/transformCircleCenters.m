function [centers,jac] = transformCircleCenters(car, x, y, th)
centers = zeros(4, 2);
jac = zeros(8, 3);

% center in map frame:
R = [cos(th), -sin(th); sin(th), cos(th)];
t = [x, y];
for i=1:4
    jac(2*(i-1)+1, :) = [1, 0, -R(2,:)*car.centers(i, :)'];
    jac(2*(i-1)+2, :) = [0, 1, R(1,:)*car.centers(i, :)'];
    centers(i, :) = car.centers(i, :) * R' + t;
end