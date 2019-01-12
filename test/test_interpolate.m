clc;
clear all;
addpath('../lib');
ctrlp = csvread('smooth_result.csv',1,0);

% the number of control points
ctrl_point_num = length(ctrlp);
% multiplicity of breakpoints
multi = 1;
% degree of spline
k = 8;
% the number of breakpoints
bpts_num = (ctrl_point_num + k - 2*k) / multi + 2;
% the number of segments
segs = bpts_num - 1;
% the configuration of specified spline
cfg = getSplineCfg(k, segs, multi);

t = 0:0.002:1.0;
% B-spline curve
B = getSplineMatrix(1, t, cfg);
x = B * ctrlp(:, 1);
y = B * ctrlp(:, 2);

dB = diffSplineMatrix(1, cfg, t);
dx = dB * ctrlp(:, 1);
dy = dB * ctrlp(:, 2);

dB2 = getDSplineMatrix(1, 1, t, cfg) * differentialMatrix(1, cfg.nc, 1);
dx2 = dB2 * ctrlp(:, 1);
dy2 = dB2 * ctrlp(:, 2);

ddB = diffSplineMatrix(2, cfg, t);
ddx = ddB * ctrlp(:, 1);
ddy = ddB * ctrlp(:, 2);

% non-uniform spline curve
w = ones(length(ctrlp), 1);
w(22:29) = [2;2;3;3;0.5;0.5;10;0.1];
R = zeros(length(t), cfg.nc);
dR = zeros(length(t), cfg.nc);
ddR = zeros(length(t), cfg.nc);
dSpline = zeros(3, cfg.nc);
for i=1:length(t)
    dSpline(1, :) = B(i, :);
    dSpline(2, :) = dB(i, :);
    dSpline(3, :) = ddB(i, :);
    param = nonUniformSplineParam(2, dSpline, w);
    R(i, :) = param(1).R;
    dR(i, :) = param(2).R;
    ddR(i, :) = param(3).R;
end
x_new = R * ctrlp(:, 1);
y_new = R * ctrlp(:, 2);
dx_new = dR * ctrlp(:, 1);
dy_new = dR * ctrlp(:, 2);
ddx_new = ddR * ctrlp(:, 1);
ddy_new = ddR * ctrlp(:, 2);


figure('Name', 'first direvative w.r.t parameters');
plot(t, dx, 'b', t, dy, 'r', t, dx_new, 'g', t, dy_new, 'y');

figure('Name', 'second direvative w.r.t parameters');
plot(t, ddx, 'b', t, ddy, 'r', t, ddx_new, 'g', t, ddy_new, 'y');

for i=1:length(t)
    k(i) = (dx(i)*ddy(i) - ddx(i)*dy(i)) / (dx(i)^2+dy(i)^2)^1.5;
    if k(i) > 10
        k(i) = 10;
    end
    k_new(i) = (dx_new(i)*ddy_new(i) - ddx_new(i)*dy_new(i)) / (dx_new(i)^2+dy_new(i)^2)^1.5;
    if k_new(i) > 10
        k_new(i) = 10;
    end
end
figure('Name', 'curvature');
plot(t, k, 'b', t, k_new, 'r');

figure('Name', 'curve')
plot(ctrlp(:,1), ctrlp(:,2), 'r.');
hold on;
plot(x, y, 'r', x_new, y_new, 'b');



