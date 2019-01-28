clc;
clear all;
addpath('../../lib');
[gausst, gaussw] = getGaussParam(50);
tau = 0.5*10*(gausst + ones(1, length(gausst)));
% N = 50;
% tau = 0:10/N:10;
%% spline-parameterized State variables
% the configuration of spline
cfg_x = getSplineCfg(6, 3, 3);
cfg_x.knot = cfg_x.knot * 10;
cfg_x.bpt = cfg_x.bpt * 10;
cfg_u = getSplineCfg(5, 4, 3);
cfg_u.knot = cfg_u.knot * 10;
cfg_u.bpt = cfg_u.bpt * 10;
% the dimension of state variable
dim_x = 3;
% the dimension of control variable
dim_u = 2;
% the number of decision variable
np = dim_x * cfg_x.nc + dim_u * cfg_u.nc;
% spline parameter matrix for state variable
B_x = getSplineMatrix(1, tau, cfg_x);
dB_x = getDSplineMatrix(1, 1, tau, cfg_x);
dP = differentialMatrix(1, cfg_x.nc, 1);
% spline parameter matrix for control variable
B_u = getSplineMatrix(1, tau, cfg_u);

%% set optimization problem
% cost function
cost = @(p)constFunc(p, B_u, cfg_x.nc, cfg_u.nc, gaussw);
% nonlinear constraints
nlon = @(p)nlonConstraints(p, dB_x*dP, B_x, B_u, cfg_x.nc, cfg_u.nc, tau);
% equality constraints
Aeq = zeros(6, np);
Aeq(1, 1) = 1;
Aeq(2, cfg_x.nc+1) = 1;
Aeq(3, 2*cfg_x.nc+1) = 1;
Aeq(4, cfg_x.nc) = 1;
Aeq(5, 2*cfg_x.nc) = 1;
Aeq(6, 3*cfg_x.nc) = 1;
beq = [10; 2; 20; 27; 8; 35];
% initial guess
p0 = zeros(np, 1);
dx = ([27;8;35]-[10;2;20]) / (cfg_x.nc-1);
for i=1:cfg_x.nc
    p0(i) = 10 + dx(1)*(i-1);
    p0(i+cfg_x.nc) = 2 + dx(2)*(i-1);
    p0(i+cfg_x.nc*2) = 20 + dx(3)*(i-1);
end
% solve
tic
snseti ('Major Iteration limit', 250);
[p,fval] = snsolve(cost, p0,[],[],Aeq,beq,[],[],nlon);
snprint off;
snend;
toc
%% plot
px1 = p(1 : cfg_x.nc);
px2 = p(cfg_x.nc+1 : 2*cfg_x.nc);
px3 = p(2*cfg_x.nc+1 : 3*cfg_x.nc);
pu1 = p(3*cfg_x.nc+1 : 3*cfg_x.nc+cfg_u.nc);
pu2 = p(3*cfg_x.nc+cfg_u.nc+1 : length(p));

figure;
plot3(B_x * px1, B_x * px2, B_x * px3);
grid on;
figure;
plot(B_u*pu1, B_u*pu2, dB_x*dP*px1, dB_x*dP*px2);


