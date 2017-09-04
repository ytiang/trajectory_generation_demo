clc;
clear all;
close all;

ref_path = getRefPath(1.0);
discrete_num = length(ref_path);
obs = obsMap();

%% Set function
dim = 3;
np = discrete_num * dim;
nlon = @(p)nlonConstraints(p, ref_path(:, 4), ref_path(:, 5), ref_path(:, 1:2), ref_path(:, 3), obs);
cost = @(p)constFunc(p, ref_path(:, 5));

%% endpoint equation constraints,Aep * x = beq
Aeq = zeros(2, np);
%d0
Aeq(1,1) = 1;
%phi0
Aeq(2,discrete_num+1) = 1;
beq = [0; -15 * pi / 180];

%% linear nonequation constraints : A*x < b
A = [];
b = [];

%% lower and upper bound:
lb = -inf*ones(np, 1);
% d low
lb(1:discrete_num) = -15*ones(discrete_num, 1);
% phi low
lb(discrete_num+1 : 2*discrete_num) = -pi / 4 * ones(discrete_num,1);
% kappa low
lb(2*discrete_num+1 : 3*discrete_num) = -0.2 * ones(discrete_num,1);
ub = inf*ones(np, 1);
ub(1:discrete_num) = 15*ones(discrete_num, 1);
% phi low
ub(discrete_num+1 : 2*discrete_num) = pi / 4 * ones(discrete_num,1);
% kappa low
ub(2*discrete_num+1 : 3*discrete_num) = 0.2 * ones(discrete_num,1);

%% initial guess
% use reference path as initial path
p0 = zeros(np, 1);

%% solve
snscreen on;
% snprint('toymin.out');  % By default, screen output is off;
sntoy.spc = which('sntoy.spc');
snspec (sntoy.spc);
snseti ('Major Iteration limit', 250);
tic
[p,fval] = snsolve(cost, p0,A,b,Aeq,beq,lb,ub,nlon);
time = toc
snprint off;
snend;


%% plot
plot_data;







