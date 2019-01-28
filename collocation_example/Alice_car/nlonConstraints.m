function [c, ceq, dc, dceq] = nlonConstraints(p, B, dB, n, T)
c = [];
dc = [];
ceq = zeros(length(T)*3, 1);
dceq = zeros(length(T)*3, length(p));
[px, py, pth, pphi] = getValue(p, n);
tf = p(length(p));
L = 2.86;
v = 1;
for i=1:length(T)
    dxi = dB.pt(i,:) * px / tf;
    dyi = dB.pt(i,:) * py / tf;
    dthi = dB.th(i,:) * pth / tf;
    thi = B.th(i,:) * pth;
    phii = B.u(i,:) * pphi;
    ceq((i-1)*3+1) = dxi - v*cos(thi);
    ceq((i-1)*3+2) = dyi - v*sin(thi);
    ceq((i-1)*3+3) = dthi - v/L*tan(phii);
    
    g1 = zeros(1, length(p));
    g1(1:n(1)) = dB.pt(i,:) / tf;
    g1(sum(n(1:2))+1 : sum(n(1:3))) = -v*(-sin(thi) * B.th(i, :));
    g1(length(p)) = -dxi / tf;
    
    g2 = zeros(1, length(p));
    g2(n(1)+1 : sum(n(1:2))) = dB.pt(i,:) / tf;
    g2(sum(n(1:2))+1 : sum(n(1:3))) = -v*(cos(thi) * B.th(i, :));
    g2(length(p)) = -dyi / tf;
    
    g3 = zeros(1, length(p));
    g3(sum(n(1:2))+1 : sum(n(1:3))) = dB.th(i,:) / tf;
    g3(sum(n(1:3))+1 : sum(n(1:4))) = -v/L*(1+(tan(phii))^2) * B.u(i, :);
    g3(length(p)) = -dthi / tf;
    
    dceq((i-1)*3+1, :) = g1;
    dceq((i-1)*3+2, :) = g2;
    dceq((i-1)*3+3, :) = g3;
end