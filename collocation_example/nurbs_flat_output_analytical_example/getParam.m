function [pz1, pz2, wz1, wz2] = getParam(p, num)
pz1 = p(1 : num);
pz2 = p(num+1 : 2*num);
wz1 = p(2*num+1 : 3*num);
wz2 = p(3*num+1 : 4*num);