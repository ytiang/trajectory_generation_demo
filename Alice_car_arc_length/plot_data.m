[px, py, pth] = getValue(p, p_num);
sf = p(np);
s_new = s(1):0.01:s(Nt);
x = SplineApprox(px, 1, s_new, cfg_pt);
y = SplineApprox(py, 1, s_new, cfg_pt);
bpt_x = SplineApprox(px, 1, s, cfg_pt);
bpt_y = SplineApprox(py, 1, s, cfg_pt);
theta = SplineApprox(pth, 1, s_new, cfg_th);

figure('Name', 'position');
plot(x, y);
hold on;
plot(px, py, '+');
plot(bpt_x, bpt_y, '*')
axis equal;

figure('Name', 'heading');
plot(s_new, theta);

coef_kappa = diffSplineMatrix(1, cfg_th, s_new);
kappa = coef_kappa * pth / sf;
figure('Name', 'curvature');
plot(s_new, kappa);