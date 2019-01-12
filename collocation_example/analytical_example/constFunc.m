function [d,g] = constFunc(p, M, x_num, u_num, weight)
pu1 = p(x_num*3+1 : x_num*3+u_num);
pu2 = p(x_num*3+u_num+1 : length(p));
d = 0;
g = zeros(length(p), 1);
n = length(weight);
w1 = 1;
w2 = 1;
for i=1:n
    u1 = M(i,:) * pu1;
    u2 = M(i,:) * pu2;
    d = d + weight(i)*(u1^2 + 2*u2^2);
    g(x_num*3+1 : x_num*3+u_num) = g(x_num*3+1 : x_num*3+u_num) + w1*weight(i)*2*u1*(M(i,:)');
    g(x_num*3+u_num+1 : length(p)) = g(x_num*3+u_num+1 : length(p)) + w2*weight(i)*2*u2*(M(i,:)');
end
d = d * 10 /4;
g = g * 10 / 4;
% g(x_num*3+1 : x_num*3+u_num) = g(x_num*3+1 : x_num*3+u_num) * 10 / 4;
% g(x_num*3+u_num+1 : length(p)) = g(x_num*3+u_num+1 : length(p)) * 10 / 4;

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