function [px, py, pth, pphi] = getValue(p, n)
px = p(1:n(1));
py = p(n(1)+1:sum(n(1:2)));
pth = p(sum(n(1:2))+1:sum(n(1:3)));
pphi = p(sum(n(1:3))+1:sum(n(1:4)));