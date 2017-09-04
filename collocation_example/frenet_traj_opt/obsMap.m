function [obs_field] = obsMap()
image = imread('obs3.png');
obs = uint8(255 * ones(size(image(:,:,1)))) - image(:,:,1);
dist = bwdist(obs);
inverse_dist = bwdist(image(:,:,1));
sdf = dist - inverse_dist;

[rows, cols] = size(sdf);

dis_threshold = 1.5;

obs_field.max_cost = 0;
obs_field.resolution = 0.1;
obs_field.rows = rows;
obs_field.cols = cols;
obs_field.cost_map = zeros(rows, cols);
obs_field.image = image;

sdf = sdf * obs_field.resolution;

obs_field.cost_map = sdf;

% for i= 1:rows
%     for j=1:cols
%         if sdf(i, j) <= 0
%             obs_field.cost_map(i, j) = dis_threshold - 2*sdf(i, j);
%         else 
%             if sdf(i, j) <= dis_threshold
%                 obs_field.cost_map(i, j) = 1 / dis_threshold * (sdf(i, j) - dis_threshold)^2;
%             end
%         end
%         obs_field.cost_map(i, j)  =min(obs_field.cost_map(i, j), 50);
%         if obs_field.cost_map(i, j) > obs_field.max_cost
%             obs_field.max_cost = obs_field.cost_map(i, j);
%         end
%     end
% end

