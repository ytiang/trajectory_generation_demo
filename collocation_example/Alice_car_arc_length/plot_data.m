[px, py, pth] = getValue(p, p_num);
sf = p(np);
s_new = s(1):0.01:s(Nt);
x = SplineApprox(px, 1, s_new, cfg_pt);
y = SplineApprox(py, 1, s_new, cfg_pt);
bpt_x = SplineApprox(px, 1, s, cfg_pt);
bpt_y = SplineApprox(py, 1, s, cfg_pt);
theta = SplineApprox(pth, 1, s_new, cfg_th);

figure('Name', 'image');
image = addPointOnImage(obs.image, bpt_x, bpt_y, obs);
imshow(image)

figure('Name', 'position');
plot(x, y);
hold on;
plot(px, py, '+');
plot(bpt_x, bpt_y, '*')
axis equal;
% check constraint:
[yc, yceq, ydc, ydceq] = nlon(p);
yobs_cost = zeros(length(s), 1);
yg_cost = zeros(length(s), 2);
for i=1:length(s)
    [a, b] = distanceField(obs, bpt_x(i), bpt_y(i));
    yobs_cost(i) = a;
    yg_cost(i,:) = b';
end

% figure('Name', 'heading');
% plot(s_new, theta);
% 
% coef_kappa = diffSplineMatrix(1, cfg_th, s_new);
% kappa = coef_kappa * pth / sf;
% figure('Name', 'curvature');
% plot(s_new, kappa);