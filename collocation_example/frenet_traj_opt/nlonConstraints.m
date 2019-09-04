function [c, ceq, dc, dceq] = nlonConstraints(p, car, ref_kappa, ref_s, ref_pt, ref_th, obs)
num = length(p) / 3;

% environment constraint:
 n = 4;
c = zeros(n*(num - 1), 1);
dc = zeros(n*(num - 1), length(p));

% vehicle kinematic constraints based on frenet:
% d(s) is the lateral offset; phi(s) = heading_traj - heading_ref
% d' = (1 - d*k_r)*tan(phi)
% phi' = k_p(1-d*k_r)/cos(phi) - k_r;
% k_r is the curvature os reference path, k_p is the curvature of actual
% path as control input.


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
    th = ref_th(i) + phi(i);
    part_d = [-sin(ref_th(i)); cos(ref_th(i)); 0];
    part_phi = [0; 0; 1];
    [centers,jac] = transformCircleCenters(car, pt(1), pt(2), th);
    for j = 1:n
        if j == 1 || j == 4
            min_r = car.r+0.3;
        else
            min_r = 0.3+car.r;
        end
        [radius, g] = distanceField(obs, centers(j,1), centers(j,2));
        c((i-2)*n + j) = min_r - radius;
        % df/d(d_i)
        dc((i-2)*n + j, i) = -g' * jac(2*j-1:2*j, :) * part_d;
        % df/d(phi_i)
        dc((i-2)*n + j, num+i) = -g' *jac(2*j-1:2*j, :) * part_phi;
    end
end