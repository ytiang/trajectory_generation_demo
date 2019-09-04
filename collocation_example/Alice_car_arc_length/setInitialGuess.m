function [p0, endpoint] = setInitialGuess(param_num)
np = sum(param_num(1:3))+1;
p0 = zeros(np, 1);
endpoint = zeros(6, 1);
curve = csvread('initial_guess.csv',1, 0);
info = zeros(length(curve), 2);
s = 0;
for i = 1:length(curve)
    if i > 1
        dp = curve(i, :) - curve(i-1, :);
        s = s + norm(dp);
    end
    if i == length(curve)
        info(i, :) = info(i-1, :);
    else
        dp = curve(i+1, :) - curve(i, :);
        info(i, 1) = atan2(dp(2), dp(1));
        info(i, 2) = s;
    end
end

endpoint(1) = curve(1,1);
endpoint(2) = curve(1,2);
endpoint(3) = info(1,1);
endpoint(4) = curve(length(curve),1);
endpoint(5) = curve(length(curve),2);
endpoint(6) = info(length(curve),1);

pt = zeros(length(param_num(1)), 2);
pt = curve(1, :);
id_pt = 2;

delta_s = s / param_num(1);

for i = 1:length(curve)
    dp = curve(i, :) - pt(id_pt-1, :);
    if(norm(dp) >= delta_s)
        pt(id_pt, :) = curve(i, :);
        id_pt = id_pt + 1;
    end
end

delta_heading = floor(length(curve) / param_num(3));
th = zeros(length(param_num(3)), 1);
for i = 1:param_num(3)
    j = (i-1) * delta_heading + 1;
    th(i) = info(j, 1);
end
th(param_num(3)) = info(length(info), 1);

p0(1:param_num(1)) = pt(:, 1);
p0(param_num(1)+1:sum(param_num(1:2))) = pt(:, 2);
p0(sum(param_num(1:2))+1:sum(param_num(1:3))) = th;
p0(np) = s;
