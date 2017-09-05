function [cost, g] = constFunc(p, B, dB, n, gaussw)
cost = 0;
g = zeros(length(p), 1);
N = length(gaussw);
id_th = sum(n(1:2)) + 1;
% pth = p(id_th:sum(n(1:3)));
[px, py, pth] = getValue(p, n);
sf = p(length(p));
gth = zeros(n(3), 1);

for i=1:N
    kappa = dB.th(i, :) * pth / sf;
    x = B.pt(i, :) * px;
    y = B.pt(i,:) * py;
    
    cost = cost + sf / 2 * gaussw(i)*(kappa^2);
    
    gth = gth + sf/2 * 2 * gaussw(i) * (kappa * (dB.th(i, :))');
    
    g(length(p)) = g(length(p)) - 1/2 * gaussw(i) * (kappa^2);
end
cost = cost + 5*sf;
g(id_th:sum(n(1:3))) = gth;
g(length(p)) = g(length(p)) + 5;
