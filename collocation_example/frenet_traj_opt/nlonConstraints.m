function [c, ceq, dc, dceq] = nlonConstraints(p, ref_kappa, ref_s, ref_pt, ref_th, obs)
num = length(p) / 3;
dis_thereshold = 1.8;

% environment constraint:
c = zeros(num - 1, 1);
dc = zeros(num - 1, length(p));
% vehicle kinematic constraints:
ceq = zeros((num - 1)*2, 1);
dceq = zeros((num - 1)*2, length(p));

d = p(1:num);
phi = p(num+1 : 2*num);
kappa = p(2*num+1 : 3*num);

for i=2:num
    ds = ref_s(i) - ref_s(i-1);
    r = (1 - d(i-1)*ref_kappa(i-1));
    derive_d = d(i-1) + ds * r * tan(phi(i-1));
    derive_phi = phi(i-1) + ds * (kappa(i-1) * r / cos(phi(i-1)) - ref_kappa(i-1) );
    
    ceq((i-2) * 2 + 1) = d(i) - derive_d;
    ceq((i-2) * 2 + 2) = phi(i) - derive_phi;
    
    dceq((i-2) * 2 + 1, i-1) = -(1 - ds*tan(phi(i-1))*ref_kappa(i-1));
    dceq((i-2) * 2 + 1, i) = 1;
    dceq((i-2) * 2 + 1, num+i-1) = -ds*r/cos(phi(i-1))^2;
    
    dceq((i-2) * 2 + 2, i-1) = -(-kappa(i-1)*ref_kappa(i-1) / cos(phi(i-1)));
    dceq((i-2) * 2 + 2, num+i-1) = -(1 + kappa(i-1)*r*tan(phi(i-1))/cos(phi(i-1)));
    dceq((i-2) * 2 + 2, num+i) = 1;
    dceq((i-2) * 2 + 2, 2*num+i-1) = -(r / cos(phi(i)));
    
    pt = ref_pt(i, :) + d(i) * [-sin(ref_th(i)), cos(ref_th(i))];
    [radius, g] = distanceField(obs, pt(1), pt(2));
    c(i-1) = dis_thereshold - radius;
    dc(i-1, i) = -[-sin(ref_th(i)), cos(ref_th(i))] * g;
end