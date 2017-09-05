function [cost, g] = distanceField(obs_map, x, y)
cost = 0;
g = zeros(2, 1);
idx = floor(x / obs_map.resolution) + 1;
idy = floor(y / obs_map.resolution) + 1;
valid = 1;
if idx < 1
    cost = cost - idx;
    g(1) = - 1.0 / obs_map.resolution;
    valid = 0;
end
if idx > obs_map.rows
    cost = cost + idx - obs_map.rows;
    g(1) = 1.0 / obs_map.resolution;
    valid = 0;
end

if idy < 1
    cost = cost - idy;
    g(2) = - 1.0 / obs_map.resolution;
    valid = 0;
end
if idy > obs_map.cols
    cost = cost + idy - obs_map.cols;
    g(2) = 1.0 / obs_map.resolution;
    valid = 0;
end

if valid == 1
    cost = obs_map.cost_map(idx, idy);
    idx_min = max(1, idx - 1);
    idx_max = min(obs_map.rows, idx + 1);
    g(1) = (obs_map.cost_map(idx_max, idy) - obs_map.cost_map(idx_min, idy)) ... 
        / (obs_map.resolution * (idx_max - idx_min)); 
    
    idy_min = max(1, idy - 1);
    idy_max = min(obs_map.cols, idy + 1);
    g(2) = (obs_map.cost_map(idx, idy_max) - obs_map.cost_map(idx, idy_min)) ... 
        / (obs_map.resolution * (idy_max - idy_min)); 
end