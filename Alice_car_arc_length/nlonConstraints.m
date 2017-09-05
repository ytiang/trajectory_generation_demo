function [c, ceq, dc, dceq] = nlonConstraints(p, B, dB, n, s)
c = zeros(length(s)*2, 1);
dc = zeros(length(s)*2, length(p));

ceq = zeros(length(s)*2, 1);
dceq = zeros(length(s)*2, length(p));
[px, py, pth] = getValue(p, n);
sf = p(length(p));

for i=1:length(s)

    dxi = dB.pt(i,:) * px / sf;
    dyi = dB.pt(i,:) * py / sf;
    thi = B.th(i,:) * pth;
    ceq((i-1)*2+1) = dxi - cos(thi);
    ceq((i-1)*2+2) = dyi - sin(thi);
    
    g1 = zeros(1, length(p));
    g1(1:n(1)) = dB.pt(i,:) / sf;
    g1(sum(n(1:2))+1 : sum(n(1:3))) = -(-sin(thi) * B.th(i, :));
    g1(length(p)) = -dxi / sf;
    
    g2 = zeros(1, length(p));
    g2(n(1)+1 : sum(n(1:2))) = dB.pt(i,:) / sf;
    g2(sum(n(1:2))+1 : sum(n(1:3))) = -(cos(thi) * B.th(i, :));
    g2(length(p)) = -dyi / sf;
    
    dceq((i-1)*2+1, :) = g1;
    dceq((i-1)*2+2, :) = g2;

    % inequality
    kappa = dB.th(i, :) * pth / sf;
    c((i-1)*2+1) = kappa - 0.2;
    c((i-1)*2+2) = -0.2 - kappa;
    
    gc1 = zeros(1, length(p));
    gc1(sum(n(1:2))+1 : sum(n(1:3))) = dB.th(i, :) / sf;
    gc1(length(p)) = - kappa / sf;
    
    gc2 = zeros(1, length(p));
    gc2(sum(n(1:2))+1 : sum(n(1:3))) = -dB.th(i, :) / sf;
    gc2(length(p)) = kappa / sf;
    
    dc((i-1)*2+1, :) = gc1;
    dc((i-1)*2+2, :) = gc2;
end