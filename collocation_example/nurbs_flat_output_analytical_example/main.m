clc;
clear all;
addpath('../../lib');
[gausst, gaussw] = getGaussParam(50);
tau = 0.5*10*(gausst + ones(1, length(gausst)));
% N = 50;
% tau = 0:10/N:10;
%% spline-parameterized State variables
% the configuration of spline
cfg = getSplineCfg(6, 4, 1);
cfg.knot = cfg.knot * 10;
cfg.bpt = cfg.bpt * 10;
% the dimension of flat output variables
dim = 4;
% the number of decision variable
np = dim * cfg.nc;
% spline parameter matrix for flat output variables
B = getSplineMatrix(1, tau, cfg);
% dB= diffSplineMatrix(1, cfg, tau);
dB = getDSplineMatrix(1, 1, tau, cfg) * differentialMatrix(1, cfg.nc, 1);
ddB = diffSplineMatrix(2, cfg, tau);

%% set optimization problem
% cost function
cost = @(p)constFunc(p, B, dB, ddB, cfg.nc, gaussw);
% nonlinear constraints
nlon = @(p)nlonConstraints(p, dB, B, cfg.nc);
% initial guess
initial_guess;
% solve
tic
% snscreen on;
% snprint('toymin.out');  % By default, screen output is off;
snseti ('Major Iteration limit', 250);
[p,fval] = snsolve(cost, p0,[],[],[],[],[],[],nlon);
snend;
snprint off;
toc
%% plot
plot_variable;



