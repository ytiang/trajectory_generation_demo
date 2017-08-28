clc;
clear all;
addpath('../../lib');
%% set problem paramters
seg_num = 4;
cfg = getSplineCfg(6, seg_num, 4);
% decesion variables degree
np = cfg.nc * 2;
% circle feasible region
region = zeros(seg_num, 3);
% circle 1: (0, 0, 3)
region(1, 1) = -1;
region(1, 2) = -1;
region(1, 3) = 3;
% circle 2: (3/2*sqrt(2), 3/2*sqrt(2), 4)
region(2, 1) = 3/2*sqrt(2);
region(2, 2) = 3/2*sqrt(2);
region(2, 3) = 4;
% circle 3: (7, 3/2*sqrt(2), 3)
region(3, 1) = 7;
region(3, 2) = 3/2*sqrt(2);
region(3, 3) = 3;
% circle 4: (11, 0, 4)
region(4, 1) = 11;
region(4, 2) = 0;
region(4, 3) = 3.5;
%% set problem
% target function:
cost = @(params)volumeCost(params, cfg);
% non-linear constraints function
nlon = @(params)nlonConstraints(params, cfg, region);
% linear equality constraints
Aeq = zeros(4, np);
Aeq(1, 1) = 1;
Aeq(2, cfg.nc+1) = 1;
Aeq(3, cfg.nc) = 1;
Aeq(4, np) = 1;
beq = [-4; -1; 14.5; 0];
% initial guess
p0 = 1:np;
p0 = p0';
%% solve inner approximation of feasile convex set
tic
snseti ('Major Iteration limit', 250);
[p,fval] = snsolve(cost, p0,[],[],Aeq,beq,[],[],nlon);
snprint off;
snend;
toc
%% solve smooth problem
% [gausst, gaussw] = getGaussParam(30);
% tau = 0.5*(gausst + ones(1, length(gausst)));
% B = getSplineMatrix(1, tau, cfg);
% cost2 = @(w)smoothCost(w, p, B, cfg);
% w0 = ones(2*cfg.nc, 1);
% tic
% snseti ('Major Iteration limit', 250);
% A = -eye(2*cfg.nc);
% b = zeros(2*cfg.nc, 1);
% [w,fval2] = snsolve(cost2, w0, A, b);
% snprint off;
% snend;
% toc
% % plot
plot_data;
