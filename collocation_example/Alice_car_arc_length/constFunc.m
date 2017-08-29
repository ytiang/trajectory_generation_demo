function [cost, g] = constFunc(p, dB, n, gaussw)
cost = 0;
g = zeros(length(p), 1);
N = length(gaussw);
id_th = sum(n(1:2))+1;
pth = p(id_th:sum(n(1:3)));
sf = p(length(p));

for i=1:N
    kappa = dB.th(i, :) * pth / sf;
    cost = cost + sf / 2 * gaussw(i)*(kappa^2);
    g(id_th:sum(n(1:3))) = g(id_th:sum(n(1:3))) + ... 
        sf/2 * 2 * gaussw(i) * (kappa * (dB.th(i, :))');
    g(length(p)) = g(length(p)) - 1/2 * gaussw(i) * (kappa^2);
end
cost = cost + 10*sf;
g(length(p)) = g(length(p)) + 10;
