function vec = nonUniformSpline(w, spline_vec)
vec = zeros(length(w), 1);
den = w' * spline_vec;
for i=1:length(w)
    vec(i) = w(i)*spline_vec(i) / den;
end