clc;
clear all;
close all;
addpath('../../lib');
[gausst, gaussw] = getGaussParam(30);
s = 0.5*(gausst + ones(1, length(gausst)));
Nt = length(s);
obs = obsMap();

%% spline-parameterized
% the order of spline basic function for state variables
cfg_pt = getSplineCfg(8, 8, 1);
cfg_th = getSplineCfg(4, 12, 1);

B.pt = getSplineMatrix(1, s, cfg_pt);
% dB.pt = getDSplineMatrix(1, 1, s, cfg_pt) * differentialMatrix(1, cfg_pt.nc, 1);
dB.pt = diffSplineMatrix(1, cfg_pt, s);
B.th = getSplineMatrix(1, s, cfg_th);
% dB.th = getDSplineMatrix(1, 1, s, cfg_th) * differentialMatrix(1, cfg_th.nc, 1);
dB.th = diffSplineMatrix(1, cfg_th, s);

% dB.u = getDSplineMatrix(1, 1, s, cfg_u) * differentialMatrix(1, cfg_u.nc, 1);


%% Set optimized problem
dim = 3;
np = cfg_pt.nc * 2 + cfg_th.nc + 1;
p_num = [cfg_pt.nc, cfg_pt.nc, cfg_th.nc];
cost = @(p)constFunc(p, B, dB, p_num, gaussw, obs);
nlon = @(p)nlonConstraints(p, B, dB, p_num, s, obs);

% set initial guess
[p0, beq] = setInitialGuess(p_num);

% endpoint equation constraints,Aep * x = beq
Aeq = zeros(dim*2, np);
%x0
Aeq(1,1:p_num(1)) = B.pt(1,:); 
%y0
Aeq(2,p_num(1)+1:sum(p_num(1:2))) = B.pt(1,:);
%th0
Aeq(3,sum(p_num(1:2))+1:sum(p_num(1:3))) = B.th(1,:);
%xn
Aeq(4,1:p_num(1)) = B.pt(Nt,:);
%yn
Aeq(5,p_num(1)+1:sum(p_num(1:2))) = B.pt(Nt,:);
%vn
Aeq(6,sum(p_num(1:2))+1:sum(p_num(1:3))) = B.th(Nt,:);

% linear nonequation constraints : A*x < b
% bounding box:
A = zeros(Nt*4, np);
b = zeros(Nt*4, 1);
% x < x_max
A(1:Nt, 1:cfg_pt.nc) = B.pt;
b(1:Nt) = (obs.rows-1) * obs.resolution * ones(Nt, 1);
% -x < 0
A(Nt+1:2*Nt, 1:cfg_pt.nc) = -B.pt;
% y < y_max
A(2*Nt+1:3*Nt, cfg_pt.nc+1:2*cfg_pt.nc) = B.pt;
b(2*Nt+1:3*Nt) = (obs.cols-1) * obs.resolution * ones(Nt, 1);
% -y < 0
A(3*Nt+1:4*Nt, cfg_pt.nc+1:2*cfg_pt.nc) = -B.pt;

% lower and upper bound:
lb = -inf*ones(np, 1);
% th low
lb(sum(p_num(1:2))+1:sum(p_num(1:3))) = -pi*ones(cfg_th.nc, 1);
% sf low
lb(np) = sqrt((beq(4)-beq(1))^2 + (beq(5)-beq(2))^2);
ub = inf*ones(np, 1);
% y upper
% ub(sum(p_num(1:1))+1:sum(p_num(1:2))) = 160;
% th upper
ub(sum(p_num(1:2))+1:sum(p_num(1:3))) = pi*ones(cfg_th.nc, 1);
% sf upper
ub(np) = 3*lb(np);

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
plot_data;








