% function [c, ceq] = nlonConstraints(p, dB, B, num)
function [c, ceq, dc, dceq] = nlonConstraints(p, dB, B, num)
[rows, cols] = size(B);
% c = zeros(rows, 1);
% dc = [rows, length(p)];
c = [];
dc = [];
ceq = zeros(6, 1);
dceq = zeros(6, length(p));
[pz1, pz2, wz1, wz2] = getParam(p, num);

spline = zeros(2, num);
for i = 1 :rows
    spline(1, :) = B(i, :);
    spline(2, :) = dB(i, :);
    param1 = nonUniformSplineParam(1, spline, wz1);
    param2 = nonUniformSplineParam(1, spline, wz2);
    z1 = param1(1).R * pz1;
    dz1 = param1(2).R * pz1;
    z2 = param2(1).R * pz2;
    dz2 = param2(2).R * pz2;
    % inequality constraint:
%     c(i) = 0.25 - dz1;
%     dc(i, 1 : num) = -param1(2).R;
%     dc(i, 2*num+1 : 3*num) = -(param1(2).W * pz1)';
   
    % end-point constraints
    if i == 1
        % ceq
        ceq(1) =  z1 - 10;
        ceq(2) =  dz2 / dz1 - 2;
        ceq(3) =  z2 - 20;

        % d(ceq1)/[d(pz1) d(wz1)]
        dceq(1, 1:num) = param1(1).R;
        dceq(1, 2*num+1 : 3*num) = (param1(1).W * pz1)';

        % d(ceq2)/[d(pz1) d(wz1) d(pz2) d(wz2)]
        dceq(2, 1:num) = -dz2 / dz1^2 * param1(2).R;
        dceq(2, 2*num+1 : 3*num) = -dz2 / dz1^2 * (param1(2).W * pz1)';
        dceq(2, num+1 : 2*num) = param2(2).R / dz1;
        dceq(2, 3*num+1 : 4*num) = (param2(2).W * pz2)' / dz1;

        % d(ceq3)/[d(pz2) d(wz2)]
        dceq(3, num+1: 2*num) = param2(1).R;
        dceq(3, 3*num+1 : 4*num) = (param2(1).W * pz1)';
    end
    if i == rows
        ceq(4) = z1 - 27;
        ceq(5) = dz2 / dz1 - 8;
        ceq(6) = z2 - 35;

        % d(ceq4)/[d(pz1) d(wz1)]
        dceq(4, 1:num) = param1(1).R;
        dceq(4, 2*num+1 : 3*num) = (param1(1).W * pz1)';

        % d(ceq5)/[d(pz1) d(wz1) d(pz2) d(wz2)]
        dceq(5, 1:num) = -dz2 / dz1^2 * param1(2).R;
        dceq(5, 2*num+1 : 3*num) = -dz2 / dz1^2 * (param1(2).W * pz1)';
        dceq(5, num+1 : 2*num) = param2(2).R / dz1;
        dceq(5, 3*num+1 : 4*num) = (param2(2).W * pz2)' / dz1;

        % d(ceq6)/[d(pz2) d(wz2)]
        dceq(6, num+1: 2*num) = param2(1).R;
        dceq(6, 3*num+1 : 4*num) = (param2(1).W * pz1)';
    end
end



