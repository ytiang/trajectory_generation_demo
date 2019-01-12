function M = nonUniformSplineMatrix(w, spline)
M = zeros(size(spline));
for i = 1: length(spline)
    vec = nonUniformSpline(w, spline(i, :)');
    M(i, :) = vec';
end