function [path] = frenetToCartesian(d, phi, ref_pt, ref_th)
num = length(d);
path = zeros(num, 3);
for i = 1 : num
    path(i, 1:2) = ref_pt(i, :) + d(i) * [-sin(ref_th(i)), cos(ref_th(i))];
    path(i, 3) = ref_th(i) + phi(i);
end