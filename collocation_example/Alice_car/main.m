clc;
clear all;
close all;
addpath('../../lib');
[gausst, gaussw] = getGaussParam(30);
tau = 0.5*(gausst + ones(1, length(gausst)));
Nt = length(tau);
path_constraint = 0;
%% spline-parameterized
% the order of spline basic function for state variables
cfg_pt = getSplineCfg(8, 8, 1);
cfg_th = getSplineCfg(6, 6, 1);
cfg_u = getSplineCfg(4, 10, 1);

B.pt = getSplineMatrix(1, tau, cfg_pt);
% dB.pt = getDSplineMatrix(1, 1, tau, cfg_pt) * differentialMatrix(1, cfg_pt.nc, 1);
dB.pt = diffSplineMatrix(1, cfg_pt, tau);
B.th = getSplineMatrix(1, tau, cfg_th);
% dB.th = getDSplineMatrix(1, 1, tau, cfg_th) * differentialMatrix(1, cfg_th.nc, 1);
dB.th = diffSplineMatrix(1, cfg_th, tau);
B.u = getSplineMatrix(1, tau, cfg_u);
% dB.u = getDSplineMatrix(1, 1, tau, cfg_u) * differentialMatrix(1, cfg_u.nc, 1);
dB.u = diffSplineMatrix(1, cfg_u, tau);
%% Set optimized problem
dim = 4;
np = cfg_pt.nc * 2 + cfg_th.nc + cfg_u.nc + 1;
p_num = [cfg_pt.nc, cfg_pt.nc, cfg_th.nc, cfg_u.nc];
cost = @(p)constFunc(p, B, dB, p_num, gaussw);
nlon = @(p)nlonConstraints(p, B, dB, p_num, tau);
% endpoint equation constraints,Aep * x = beq
Aeq = zeros(dim*2, np);
%x0
Aeq(1,1:p_num(1)) = B.pt(1,:); 
%y0
Aeq(2,p_num(1)+1:sum(p_num(1:2))) = B.pt(1,:);
%th0
Aeq(3,sum(p_num(1:2))+1:sum(p_num(1:3))) = B.th(1,:);
%phi0
Aeq(4,sum(p_num(1:3))+1:sum(p_num(1:4))) = B.u(1,:);
%xn
Aeq(5,1:p_num(1)) = B.pt(Nt,:);
%yn
Aeq(6,p_num(1)+1:sum(p_num(1:2))) = B.pt(Nt,:);
%vn
Aeq(7,sum(p_num(1:2))+1:sum(p_num(1:3))) = B.th(Nt,:);
%thn
Aeq(8,sum(p_num(1:3))+1:sum(p_num(1:4))) = B.u(Nt,:);
beq = [0;0;0/180*pi;0; 25;25;-90/180*pi;0];

% linear nonequation constraints : A*x < b
if(path_constraint)
    A = zeros(2*cfg_pt.nc, np);
    b = 30*ones(2*cfg_pt.nc, 1);
    A(1:cfg_pt.nc, 1:2*cfg_pt.nc) = [-eye(cfg_pt.nc), eye(cfg_pt.nc)];
    A(cfg_pt.nc+1:2*cfg_pt.nc, 1:2*cfg_pt.nc) = [eye(cfg_pt.nc), -eye(cfg_pt.nc)];
else
    A = [];
    b = [];
end

% lower and upper bound:
lb = -inf*ones(np, 1);
% th low
lb(sum(p_num(1:2))+1:sum(p_num(1:3))) = -pi*ones(cfg_th.nc, 1);
% phi low
lb(sum(p_num(1:3))+1:sum(p_num(1:4))) = -25/180*pi*ones(cfg_u.nc,1);
% tf low
lb(np) = 0.1;
ub = inf*ones(np, 1);
% y upper
% ub(sum(p_num(1:1))+1:sum(p_num(1:2))) = 160;
% th upper
ub(sum(p_num(1:2))+1:sum(p_num(1:3))) = pi*ones(cfg_th.nc, 1);
% phi upper
ub(sum(p_num(1:3))+1:sum(p_num(1:4))) = 25/180*pi*ones(cfg_u.nc,1);
% tf upper
ub(np) = 200;
% initial guess:
p0 = zeros(np, 1);
dx = (beq(5:6)-beq(1:2))/(cfg_pt.nc-1);
for i=1:cfg_pt.nc
    p0(i) = (i-1)*dx(1);
    p0(i+cfg_pt.nc) = (i-1)*dx(2);
end
dth = (beq(7)-beq(3))/(cfg_th.nc-1);
for i=1:cfg_th.nc
    p0(i+sum(p_num(1:2))) = (i-1)*dth;
end
p0(np) = 10;
% options = optimoptions(@fmincon,'Algorithm','sqp-legacy','MaxIterations',3000);
% [p,fval] = fmincon(cost, p0,A,b,Aeq,beq,lb,ub,nlon, options);

% snscreen on;
% snprint('toymin.out');  % By default, screen output is off;
sntoy.spc = which('sntoy.spc');
snspec (sntoy.spc);
snseti ('Major Iteration limit', 250);
tic
[p,fval] = snsolve(cost, p0,A,b,Aeq,beq,lb,ub,nlon);
time = toc
% snprint off;
snend;
%% plot
[px, py, pth, pphi] = getValue(p, p_num);
tf = p(np);
t = tau(1):0.01:tau(Nt);
x = SplineApprox(px, 1, t, cfg_pt);
y = SplineApprox(py, 1, t, cfg_pt);
pth = SplineApprox(pth, 1, t, cfg_th);
phi = SplineApprox(pphi, 1, t, cfg_u);
figure('Name', 'position');
plot(x, y);
hold on;
plot(px, py, '+');
bpt_x = SplineApprox(px, 1, tau, cfg_pt);
bpt_y = SplineApprox(py, 1, tau, cfg_pt);
plot(bpt_x, bpt_y, '*')
axis equal;
if(path_constraint)
    bdx = 0:10:150;
    plot(bdx, bdx+30*ones(1, length(bdx)), bdx, bdx-30*ones(1, length(bdx)));
end
figure('Name', 'heading');
plot(tf*t, pth*180/pi);
figure('Name', 'steering');
plot(tf*t, phi*180/pi);




