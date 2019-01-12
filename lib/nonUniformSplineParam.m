function param = nonUniformSplineParam(order, spline, w)
[rows, cols] = size(spline);
if rows < order + 1
    error('order %d NURB requires input %d_th derivative of B spline vector', order, order+1);
end
% param(i) member: (i means the i-th derivative)
%     den  : w' * [B_i(t)]
%     diag : diag{[B_i(t)]}; a diagonal matrix constructed from B-spline vector
%     R    : the NURB base vector, which is a row vector
%     W    : direvative matrix of the NURB vector w.r.t weigth vector w;
cfg = struct('den', 0, ...
    'diag', zeros(length(w)), ...
    'R', zeros(1, length(w)), ...
    'W', zeros(length(w)));
param = repmat(cfg, order + 1, 1);
F0 = spline(1, :)' * w' / (spline(1, :) * w) - eye(length(w));
for i = 1 : order+1
%     w
%     bi = spline(i, :)
%     yt_debug = bi * w
    param(i).den = spline(i, :) * w;
    param(i).diag = diag(spline(i, :));
    param(i).R = w' * param(i).diag / param(1).den;
    param(i).W = -F0 * param(i).diag / param(1).den;
    r = i - 1;
    for k = 1:r
        c = nchoosek(r, k);
        R = c * param(k+1).den / param(1).den * param(r-k+1).R;
        W1 = c * F0 * spline(k+1, :)' * param(r-k+1).R / param(1).den;
        W2 = c * param(k+1).den / param(1).den * param(r-k+1).W;
        param(i).R = param(i).R - R;
        param(i).W = param(i).W + W1 - W2;
    end
end