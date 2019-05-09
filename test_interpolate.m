clc;
clear all;
addpath('../lib');
ctrlp = csvread('smooth_result.csv',1,0);

% the number of control points
ctrl_point_num = length(ctrlp);
% multiplicity of breakpoints
multi = 1;
% degree of spline
k = 6;
% the number of breakpoints
bpts_num = (ctrl_point_num + k - 2*k) / multi + 2;
% the number of segments
segs = bpts_num - 1;
% the configuration of specified spline
cfg = getSplineCfg(k, segs, multi);

% parameter
t = 0:0.001:1.0;
B = getSplineMatrix(1, t, cfg);
x = B * ctrlp(:, 1);
y = B * ctrlp(:, 2);

plot(ctrlp(:,1), ctrlp(:,2), 'r.');
hold on;
plot(x, y);

[dctrlp, dcfg] = differentialCtrlp(1, ctrlp, cfg);
[ddctrlp, ddcfg] = differentialCtrlp(1, dctrlp, dcfg);

dB = getSplineMatrix(1, t, dcfg);
ddB = getSplineMatrix(1, t, ddcfg);

dx = dB * dctrlp(:, 1);
dy = dB * dctrlp(:, 2);
figure
plot(t, dx, t, dy);

figure
ddx = ddB * ddctrlp(:, 1);
ddy = ddB * ddctrlp(:, 2);
plot(t, ddx, t, ddy);

for i=1:length(t)
    k(i) = (dx(i)*ddy(i) - ddx(i)*dy(i)) / (dx(i)^2+dy(i)^2)^1.5;
end
figure
plot(t, k);



