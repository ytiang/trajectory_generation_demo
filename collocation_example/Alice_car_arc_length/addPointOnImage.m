function image=addPointOnImage(src, x, y, obs)
image = src;
for i = 1:obs.rows
    for j = 1:obs.cols
        rho = 1.0 - obs.cost_map(i, j) / obs.max_cost;
        image(i, j, 1) = 255 * rho;
        image(i, j, 2) = 255 * rho;
        image(i, j, 3) = 255 * rho;
    end
end
for i = 1:length(x)
    idx = floor(x(i) / obs.resolution) + 1;
    idy = floor(y(i) / obs.resolution) + 1;
    x_start = max(1, idx - 2);
    x_end = min(obs.rows, idx+2);
    y_start = max(1, idy - 2);
    y_end = min(obs.cols, idy+2);
    for i=x_start:x_end
        for j=y_start:y_end
            image(i, j, 1) = 100;
            image(i, j, 2) = 100;
        end
    end
end
    
        