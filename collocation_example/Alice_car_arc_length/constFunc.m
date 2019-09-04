function [cost, g] = constFunc(p, B, dB, n, gaussw, obs)
cost = 0;
g = zeros(length(p), 1);
N = length(gaussw);
id_th = sum(n(1:2)) + 1;
% pth = p(id_th:sum(n(1:3)));
[px, py, pth] = getValue(p, n);
sf = p(length(p));
gth = zeros(n(3), 1);
gx = zeros(n(1), 1);
gy = zeros(n(2), 1);
w = 0;

for i=1:N
    kappa = dB.th(i, :) * pth / sf;
    x = B.pt(i, :) * px;
    y = B.pt(i,:) * py;
    [obs_cost, obs_grad] = distanceField(obs, x, y);
    
    cost = cost + sf / 2 * gaussw(i)*(kappa^2 + w*obs_cost);
    
    gth = gth + sf/2 * 2 * gaussw(i) * (kappa * (dB.th(i, :))');
    
    gx = gx + sf / 2 * 2 * gaussw(i) * (w*obs_grad(1) * B.pt(i, :)');
    gy = gy + sf / 2 * 2 * gaussw(i) * (w*obs_grad(2) * B.pt(i, :)');
    
    g(length(p)) = g(length(p)) - 1/2 * gaussw(i) * (kappa^2);
end
cost = cost + 5*sf;
g(id_th:sum(n(1:3))) = gth;
g(1: n(1)) = gx;
g(n(1) + 1 : sum(n(1:2))) = gy;
g(length(p)) = g(length(p)) + 5;
