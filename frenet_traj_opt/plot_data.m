
d = p(1:discrete_num);
phi = p(discrete_num+1 : 2*discrete_num);
kappa = p(2*discrete_num+1 : 3*discrete_num);
path = frenetToCartesian(d, phi, ref_path(:, 1:2), ref_path(:, 3));

% figure('Name', 'image');
% image = addPointOnImage(obs.image, bpt_x, bpt_y, obs);
% imshow(image)
colordef black
figure('Name', 'position');
[B,L] = bwboundaries(obs.image(:,:, 1),'noholes');
boundaries = B{1, 1}*obs.resolution;
patch(boundaries(:, 1), boundaries(:, 2), 'w');
hold on;
plot(path(:,1), path(:,2), 'r', ref_path(:,1), ref_path(:, 2), 'b');
hold on;

for i=1:discrete_num
    R = [cos(path(i, 3)), -sin(path(i, 3)); sin(path(i, 3)), cos(path(i, 3))];
    t = path(i, 1:2)';
    foot = zeros(2, 5);
    for j = 1:4
        foot(:, j) = R * car.footprint(j, :)' + t;
    end
    foot(:, 5) = foot(:, 1);
    plot(foot(1, :), foot(2, :));
    hold on;
end

% 
% figure('Name', 'cost map');
% x_range = 1:obs.rows;
% y_range = 1:obs.cols;
% [X, Y] = meshgrid(x_range, y_range);
% surf(x_range, y_range, obs.cost_map);

% figure('Name', 'heading');
% plot(s_new, theta);
% 
figure('Name', 'curvature');
plot(ref_path(:, 5), kappa);