function [c, ceq, dc, dceq] = nlonConstraints(p, dMX, MX, MU, x_num, u_num, T)
c = [];
dc = [];
ceq = zeros(length(T)*3, 1);
dceq = zeros(length(T)*3, length(p));

px1 = p(1 : x_num);
px2 = p(x_num+1 : 2*x_num);
px3 = p(2*x_num+1 : 3*x_num);
pu1 = p(3*x_num+1 : 3*x_num+u_num);
pu2 = p(3*x_num+u_num+1 : length(p));

for i=1:length(T)
    dx1 = dMX(i,:) * px1;
    dx2 = dMX(i,:) * px2;
    dx3 = dMX(i,:) * px3;
    x2 = MX(i,:) * px2;
    u1 = MU(i,:) * pu1;
    u2 = MU(i,:) * pu2;
    ceq((i-1)*3 + 1) = dx1 - u1;  
    ceq((i-1)*3 + 2) = dx2 - u2;
    ceq((i-1)*3 + 3) = dx3 - x2*u1;
    
    g1 = zeros(1, length(p));   
    g1(1:x_num) = dMX(i, :);
    g1(3*x_num+1 : 3*x_num+u_num) = -MU(i, :);
    
    g2 = zeros(1, length(p));   
    g2(x_num+1 : 2*x_num) = dMX(i, :);
    g2(3*x_num+u_num+1 : length(p)) = -MU(i, :);
    
    g3 = zeros(1, length(p));   
    g3(2*x_num+1 : 3*x_num) = dMX(i, :);
    g3(x_num+1 : 2*x_num) = -u1*MX(i, :);
    g3(3*x_num+1 : 3*x_num+u_num) = -x2 * MU(i, :);
    
    dceq((i-1)*3 + 1, :) = g1; 
    dceq((i-1)*3 + 2, :) = g2;
    dceq((i-1)*3 + 3, :) = g3;
end
