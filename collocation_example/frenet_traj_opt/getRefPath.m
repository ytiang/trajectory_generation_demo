function [ref_path] = getRefPath(delta_s)
% delta_s = 1.0;
curve = csvread('reference_path.csv',1, 0);
num = length(curve);
ref_path = [];
actual_num = 1;
for i = 1:num
    if i == 1
        ref_path(actual_num, :) = curve(i, :);
        actual_num = actual_num + 1;
    else
        ds = curve(i, 5) - ref_path(actual_num-1, 5);
        if ds >= delta_s
            ref_path(actual_num, :) = curve(i, :);
            actual_num = actual_num + 1;
        end
    end
end
