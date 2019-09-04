function [cost, g] = constFunc(p, ref_s)
num = length(p) / 3;
cost = 0;
g = zeros(length(p), 1);
kappa = p(num*2+1 : num*3);
d = p(1:num);
g_kappa = zeros(num, 1);
g_d = zeros(num, 1);
w1 = 1;
w2 = 0;
for i=1:num    
    cost = cost + kappa(i)^2;
    g_kappa(i) = 2 * kappa(i);
    if i > 1
        ds = ref_s(i) - ref_s(i-1);
        d_d = (d(i) - d(i-1)) / ds;
        d_kappa = (kappa(i) - kappa(i-1)) / ds;
        cost = cost + w1*(d_d)^2 + w2 * d_kappa^2;
        g_d(i-1) = g_d(i - 1) - 2*w1*d_d; 
        g_d(i) = g_d(i) + 2*w1*d_d; 
        g_kappa(i-1) = g_kappa(i-1) - 2*w2*d_kappa;
        g_kappa(i) = g_kappa(i) + 2*w2*d_kappa;
    end
end
g(num*2+1 : num*3) = g_kappa;
g(1:num) = g_d;
