p0 = ones(np, 1);
num = length(tau);
x1_0 = 1.7 * tau' + 10 * ones(num, 1);
u1_0 = 1.7 * ones(num, 1);
x2_0 = 0.2471 * (tau').^2 - 1.8706 * tau' + 2 * ones(num, 1);
u2_0 = 2 * 0.2471 * tau' - 1.8706;
x3_0 = 1.7 * (0.2471 / 3 * (tau').^3 - 1.8706/2 * (tau').^2 + 2 * tau') + 20;
p0(1:cfg_x.nc) = B.x \ x1_0;
p0(cfg_x.nc+1: 2*cfg_x.nc) = B.x \ x2_0;
p0(2*cfg_x.nc+1: 3*cfg_x.nc) = B.x \ x3_0;
p0(3*cfg_x.nc+1 : 3*cfg_x.nc+cfg_u.nc) = B.u \ u1_0;
p0(3*cfg_x.nc+cfg_u.nc+1 : 3*cfg_x.nc+2*cfg_u.nc) = B.u \ u2_0;

