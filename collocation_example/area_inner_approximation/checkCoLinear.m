function flag = checkCoLinear(x, y)
num = length(x);
flag =1;
index = 1;
for i = 2:num
    index = i;
    base_vec = [x(i) - x(1); y(i) - y(1)];
    if abs(base_vec(1)) > 1e-4 || abs(base_vec(2)) > 1e-4
        break;
    end
end
if index + 3 > num
    flag = 1;
    return;
end
count = 1;
for i = index+1 : num
    vec = [x(i) - x(1); y(i) - y(1)];
    if abs(vec(1)) < 1e-4 && abs(vec(2)) < 1e-4
        continue;
    end
    if abs(base_vec(1)) < 1e-4
        if abs(vec(1)) > 1e-4
            count = count + 1;
        end
    else
        k = vec(1) / base_vec(1);
        if abs(base_vec(2) * k - vec(1)) > 1e-4
            count = count + 1;
        end
    end
    if count >= 3
        flag = 0;
        return;
    end
end
flag = 1;