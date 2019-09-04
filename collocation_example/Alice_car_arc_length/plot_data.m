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
figure('Name', 'cost map');
x_range = 1:obs.rows;
y_range = 1:obs.cols;
[X, Y] = meshgrid(x_range, y_range);
surf(x_range, y_range, obs.cost_map);

% figure('Name', 'heading');
% plot(s_new, theta);
% 
% coef_kappa = diffSplineMatrix(1, cfg_th, s_new);
% kappa = coef_kappa * pth / sf;
% figure('Name', 'curvature');
% plot(s_new, kappa);