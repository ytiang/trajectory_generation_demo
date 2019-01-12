function [c, ceq, dc, dceq] = nlonConstraints(p, dB, B, x_num, u_num)
% function [c, ceq] = nlonConstraints(p, dB, B, x_num, u_num)
[rows, cols] = size(B.x);
c = [];
dc = [];
ceq = zeros(rows*3 + 6, 1);
dceq = zeros(rows*3 + 6, length(p));

[px1, px2, px3, pu1, pu2, wx1, wx2, wx3, wu1, wu2] = getParam(p, x_num, u_num);
dBx = zeros(2, x_num);
for i=1:rows
    dBx(1, :) = B.x(i, :);
    dBx(2, :) = dB.x(i, :);
    cfg_u1 = nonUniformSplineParam(0, B.u(i, :), wu1);
    cfg_u2 = nonUniformSplineParam(0, B.u(i, :), wu2);
    cfg_x1 = nonUniformSplineParam(1, dBx, wx1);
    cfg_x2 = nonUniformSplineParam(1, dBx, wx2);
    cfg_x3 = nonUniformSplineParam(1, dBx, wx3);
    dx1 = cfg_x1(2).R * px1;
    dx2 = cfg_x2(2).R * px2;
    dx3 = cfg_x3(2).R * px3;
    u1 = cfg_u1(1).R * pu1;
    u2 = cfg_u2(1).R * pu2;
    x2 = cfg_x2(1).R * px2;
    %% differential constraints
    ceq((i-1)*3 + 1) = dx1 - u1;  
    ceq((i-1)*3 + 2) = dx2 - u2;
    ceq((i-1)*3 + 3) = dx3 - x2*u1;
    
    % d(ceq1)/[d(px1), d(pu1), d(wx1), d(wu1)]
    % d(px1)
    dceq((i-1)*3 + 1, 1:x_num) = cfg_x1(2).R;
    % d(pu1)
    dceq((i-1)*3 + 1, 3*x_num+1 : 3*x_num+u_num) = -cfg_u1(1).R;
    % d(wx1)
    dceq((i-1)*3 + 1, 3*x_num+2*u_num+1 : 4*x_num+2*u_num) = (cfg_x1(2).W * px1)';
    % d(wu1)
    dceq((i-1)*3 + 1, 6*x_num+2*u_num+1 : 6*x_num+3*u_num) = -(cfg_u1(1).W * pu1)';
    
    % d(ceq2)/[d(px2), d(pu2), d(wx2), d(wu2)]
    % d(px2)
    dceq((i-1)*3 + 2, x_num+1 : 2*x_num) = cfg_x2(2).R;
    % d(pu2)
    dceq((i-1)*3 + 2, 3*x_num+u_num+1 : 3*x_num+2*u_num) = -cfg_u2(1).R;
    % d(wx2)
    dceq((i-1)*3 + 2, 4*x_num+2*u_num+1 : 5*x_num+2*u_num) = (cfg_x2(2).W * px2)';
    % d(wu2)
    dceq((i-1)*3 + 2, 6*x_num+3*u_num+1 : 6*x_num+4*u_num) = -(cfg_u2(1).W * pu2)';

    % d(ceq3)/[d(px3), d(px2), d(pu1), d(wx3), d(wx2), d(wu1)]
    % d(px3)
    dceq((i-1)*3 + 3, 2*x_num+1 : 3*x_num) = cfg_x3(2).R;
    % d(px2)
    dceq((i-1)*3 + 3, x_num+1 : 2*x_num) = -u1 * cfg_x2(1).R;
    % d(pu1)
    dceq((i-1)*3 + 3, 3*x_num+1 : 3*x_num+u_num) = -x2 * cfg_u1(1).R;
    % d(wx3)
    dceq((i-1)*3 + 3, 5*x_num+2*u_num+1 : 6*x_num+2*u_num) = (cfg_x3(2).W * px3)';
    % d(wx2)
    dceq((i-1)*3 + 3, 4*x_num+2*u_num+1 : 5*x_num+2*u_num) = -u1 * (cfg_x2(1).W * px2)';
    % d(wu1)
    dceq((i-1)*3 + 3, 6*x_num+2*u_num+1 : 6*x_num+3*u_num) = -x2 * (cfg_u1(1).W * pu1)';
    %% endpoint constraints
    if i == 1
        x1 = cfg_x1(1).R * px1;
        x3 = cfg_x3(1).R * px3;
        ceq(rows*3 + 1) = x1 - 10;
        ceq(rows*3 + 2) = x2 - 2;
        ceq(rows*3 + 3) = x3 - 20;
        
        dceq(rows*3 + 1, 1:x_num) = cfg_x1(1).R;
        dceq(rows*3 + 1, 3*x_num+2*u_num+1 : 4*x_num+2*u_num) = (cfg_x1(1).W*px1)';
        dceq(rows*3 + 2, x_num+1 : 2*x_num) = cfg_x2(1).R;
        dceq(rows*3 + 2, 4*x_num+2*u_num+1 : 5*x_num+2*u_num) = (cfg_x2(1).W*px2)';
        dceq(rows*3 + 3, 2*x_num+1 : 3*x_num) = cfg_x3(1).R;
        dceq(rows*3 + 3, 5*x_num+2*u_num+1 : 6*x_num+2*u_num) = (cfg_x3(1).W*px3)';
    end
    if i == rows
        x1 = cfg_x1(1).R * px1;
        x3 = cfg_x3(1).R * px3;
        ceq(rows*3 + 4) = x1 - 27;
        ceq(rows*3 + 5) = x2 - 8;
        ceq(rows*3 + 6) = x3 - 35;
        
        dceq(rows*3 + 4, 1:x_num) = cfg_x1(1).R;
        dceq(rows*3 + 4, 3*x_num+2*u_num+1 : 4*x_num+2*u_num) = (cfg_x1(1).W*px1)';
        dceq(rows*3 + 5, x_num+1 : 2*x_num) = cfg_x2(1).R;
        dceq(rows*3 + 5, 4*x_num+2*u_num+1 : 5*x_num+2*u_num) = (cfg_x2(1).W*px2)';
        dceq(rows*3 + 6, 2*x_num+1 : 3*x_num) = cfg_x3(1).R;
        dceq(rows*3 + 6, 5*x_num+2*u_num+1 : 6*x_num+2*u_num) = (cfg_x3(1).W*px3)';
    end
end

