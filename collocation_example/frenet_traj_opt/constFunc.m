function [cost, g] = constFunc(p, ref_s)
num = length(p) / 3;
cost = 0;
g = zeros(length(p), 1);
kappa = p(num*2+1 : num*3);
d = p(1:num);
g_kappa = zeros(num, 1);
g_d = zeros(num, 1);
w = 5;
for i=1:num    
    cost = cost + kappa(i)^2;
    g_kappa(i) = 2 * kappa(i);
    if i > 1
        ds = ref_s(i) - ref_s(i-1);
        d_d = (d(i) - d(i-1)) / ds;
        cost = cost + w*(d_d)^2;
        g_d(i-1) = g_d(i - 1) - 2*w*d_d; 
        g_d(i) = g_d(i) + 2*w*d_d; 
    end
end
g(num*2+1 : num*3) = g_kappa;
g(1:num) = g_d;
