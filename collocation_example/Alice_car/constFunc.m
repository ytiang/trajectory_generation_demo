function [cost, g] = constFunc(p, B, dB, n, gaussw)
cost = 0;
g = zeros(length(p), 1);
N = length(gaussw);
id_u = sum(n(1:3))+1;
pphi = p(id_u:sum(n(1:4)));
tf = p(length(p));
w = 0;
for i=1:N
    phi = B.u(i, :) * pphi;
    dphi = dB.u(i, :) * pphi / tf;
    cost = cost + tf / 2 * gaussw(i)*(phi^2 + w*dphi^2);
    g(id_u:sum(n(1:4))) = g(id_u:sum(n(1:4))) + ... 
        tf/2 * 2 * gaussw(i) * (phi * (B.u(i, :))' + w*dphi/tf * (dB.u(i, :))');
    g(length(p)) = g(length(p)) + 1/2 * gaussw(i) * (phi^2 - w*dphi^2 / tf / tf);
end
cost = cost + 0.01*tf;
g(length(p)) = g(length(p)) + 0.01;
