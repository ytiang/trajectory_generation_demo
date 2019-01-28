clc;
clear all;
addpath('../lib');
%%
real = 2 / pi + 1/3*(5^3-1^3) - 3*(5^2-1^2) + 5*(5-1);

%% Gauss Quadrature:
[gausst, gaussw] = getGaussParam(30);
d = 0;
x = 5 / 2 * (gausst + ones(1, length(gausst)));
for i=1:length(gaussw)
    d = d + myFunc(x(i))*gaussw(i);
end
d = d * 5 / 2;

%% simpson
N = 15;
n = 3;
dt = 5 / 30;
d2 = 0;
for i = 1:6
    x0 = (n * (i-1)) * dt;
    x1 = x0 + dt;
    x2 = x1 + dt;
    x3 = x2 + dt;
    x4 = x3 + dt;
    x5 = x4 + dt;
%     d2 = d2 + 3 / 8 * dt * (myFunc(x0) + 3 * myFunc(x1) + 3 * myFunc(x2) + myFunc(x3));
d2 = d2 + 5 / 288 * dt * (19*myFunc(x0) + 75 * myFunc(x1) + 50 * myFunc(x2) + 50*myFunc(x3) + 75*myFunc(x4) + 19*myFunc(x5));
end