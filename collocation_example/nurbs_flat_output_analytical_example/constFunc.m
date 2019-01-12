% function [d] = constFunc(p, B, dB, ddB, num, weight)
function [d,g] = constFunc(p, B, dB, ddB, num, weight)
[pz1, pz2, wz1, wz2] = getParam(p, num);
d = 0;
g = zeros(length(p), 1);
n = length(weight);
w1 = 1;
w2 = 1;
spline = zeros(3, num);
for i=1:n
%     count = count + 1
    spline(1, :) = B(i, :);
    spline(2, :) = dB(i, :);
    spline(3, :) = ddB(i, :);
    param1 = nonUniformSplineParam(2, spline, wz1);
    param2 = nonUniformSplineParam(2, spline, wz2);
    dz1 = param1(2).R * pz1;
    ddz1 = param1(3).R * pz1;
    dz2 = param2(2).R * pz2;
    ddz2 = param2(3).R * pz2;
    u1 = dz1;
    u2 = (ddz2*dz1 - ddz1*dz2) / dz1^2;
    d = d + weight(i)*(w1*u1^2 + w2*u2^2);
    % df/dpz1
    g(1 : num) = g(1 : num) + weight(i) * ... 
        ( w1*2*u1 * param1(2).R' + ...
          w2*2*u2 * ( (ddz2 * param1(2).R' - dz2 * param1(3).R') / dz1^2 - ... 
                      2* u2 / dz1 * param1(2).R' ) );
    % df/dpz2
    g(num+1 : 2*num) = g(num+1 : 2*num) + w2*weight(i)*2*u2 * ... 
        (dz1 * param2(3).R' - ddz1 * param2(2).R') / dz1^2; 
    % df/dwz1
    g(2*num+1 : 3*num) = g(2*num+1 : 3*num) + weight(i)* ... 
        ( w1*2*u1 * param1(2).W*pz1 + ... 
          w2*2*u2 * ( (ddz2*param1(2).W*pz1 - dz2 * param1(3).W*pz1) / dz1^2 - ...
                      2* u2 / dz1 * param1(2).W*pz1) );
    % df/dwz2
    g(3*num+1 : 4*num) = g(3*num+1 : 4*num) - w2 * weight(i)*2*u2 * ... 
        (dz1*param2(3).W*pz2 - ddz1*param2(2).W*pz2) / dz1^2;
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