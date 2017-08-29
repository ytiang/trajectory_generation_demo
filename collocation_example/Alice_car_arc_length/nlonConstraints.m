function [c, ceq, dc, dceq] = nlonConstraints(p, B, dB, n, s)
c = [];
dc = [];
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
end