function [d,g] = constFunc(p, B, x_num, u_num, weight)
% function [d] = constFunc(p, B, x_num, u_num, weight)
[px1, px2, px3, pu1, pu2, wx1, wx2, wx3, wu1, wu2] = getParam(p, x_num, u_num);
d = 0;
g = zeros(length(p), 1);
n = length(weight);
w1 = 1;
w2 = 1;
for i=1:n
    param1 = nonUniformSplineParam(0, B(i, :), wu1);
    param2 = nonUniformSplineParam(0, B(i, :), wu2);
    u1 = param1(1).R * pu1;
    u2 = param2(1).R * pu2;
    d = d + weight(i)*(w1 * u1^2 + w2 * u2^2);
    
    % d(pu1)
    g(3*x_num+1 : 3*x_num+u_num) = ... 
        g(3*x_num+1 : 3*x_num+u_num) + ... 
        weight(i)*w1*2*u1 * param1(1).R';
    % d(pu2)
    g(3*x_num+u_num+1 : 3*x_num+2*u_num) = ... 
        g(3*x_num+u_num+1 : 3*x_num+2*u_num) + ... 
        weight(i)*w2*2*u2*param2(1).R';
    % d(wu1)
    g(6*x_num+2*u_num+1 : 6*x_num+3*u_num) = ... 
        g(6*x_num+2*u_num+1 : 6*x_num+3*u_num) + ... 
        weight(i)*w1*2*u1*param1(1).W*pu1;
    % d(wu2)
    g(6*x_num+3*u_num+1 : 6*x_num+4*u_num) = ... 
        g(6*x_num+3*u_num+1 : 6*x_num+4*u_num) + ... 
        weight(i)*w2*2*u2*param2(1).W*pu2;
end
d = d * 10 /4;
g = g * 10 / 4;

% n = 5;
% w = 5 / 288 * 10 / 50 * [19, 75, 50, 50, 75, 19];
% for i = 1:10
%     j = (i-1)*n+1;
%     u1 = M(j:j+n, :) * pu1;
%     u2 = M(j:j+n, :) * pu2;
%     d = d + w * (u1.^2 + 2*u2.^2);
%     for l=1:n+1
%         g(x_num*3+1 : x_num*3+u_num) = g(x_num*3+1 : x_num*3+u_num) + 2*w(l)*u1(l)*(M(j+l-1, :)');
%         g(x_num*3+u_num+1 : length(p)) = g(x_num*3+u_num+1 : length(p)) + 2*2*w(l)*u2(l)*(M(j+l-1,:)');
%     end
% end