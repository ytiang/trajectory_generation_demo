[px1, px2, px3, pu1, pu2, wx1, wx2, wx3, wu1, wu2] = getParam(p, cfg_x.nc, cfg_u.nc);
x1 = zeros(length(tau), 1);
x2 = zeros(length(tau), 1);
x3 = zeros(length(tau), 1);
u1 = zeros(length(tau), 1);
u2 = zeros(length(tau), 1);
dx1 = zeros(length(tau), 1);
dx2 = zeros(length(tau), 1);
dx3 = zeros(length(tau), 1);
inter_x1 = zeros(length(tau), 1);
inter_x2 = zeros(length(tau), 1);
inter_x3 = zeros(length(tau), 1);
for i=1:length(tau)
    dBx(1, :) = B.x(i, :);
    dBx(2, :) = dB.x(i, :);
    cfg_u1 = nonUniformSplineParam(0, B.u(i, :), wu1);
    cfg_u2 = nonUniformSplineParam(0, B.u(i, :), wu2);
    cfg_x1 = nonUniformSplineParam(1, dBx, wx1);
    cfg_x2 = nonUniformSplineParam(1, dBx, wx2);
    cfg_x3 = nonUniformSplineParam(1, dBx, wx3);
    x1(i) = cfg_x1(1).R * px1;
    x2(i) = cfg_x2(1).R * px2;
    x3(i) = cfg_x3(1).R * px3;
    u1(i) = cfg_u1(1).R * pu1;
    u2(i) = cfg_u2(1).R * pu2;
    dx1(i) = cfg_x1(2).R * px1;
    dx2(i) = cfg_x2(2).R * px2;
    dx3(i) = cfg_x3(2).R * px3;
    
    if i == length(tau) 
        inter_x1(i) = 27;
        inter_x2(i) = 8;
        inter_x3(i) = 35;
    else
        if i == 1
        inter_x1(i) = 10;
        inter_x2(i) = 2;
        inter_x3(i) = 20;
        else
            dt = tau(i) - tau(i-1);
            inter_x1(i) = inter_x1(i-1) + dt*(u1(i-1)+u1(i))/2;
            inter_x2(i) = inter_x2(i-1) + dt*(u2(i-1)+u2(i))/2;
            inter_x3(i) = inter_x3(i-1) + dt*(inter_x2(i-1)*u1(i-1)+inter_x2(i)*u1(i))/2;
        end
    end
end

figure('Name', 'state')
plot3(x1, x2, x3);
grid on;
figure('Name', 'control')
plot(u1, u2);
grid on;

figure('Name', 'x');
subplot(3,1,1);
plot(tau, x1, 'b', tau, x1_0, 'k-.', tau, inter_x1, 'm.');
axis([0, 10, 10, 30]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[10:5:30]);
grid on;
subplot(3,1,2);
plot(tau, x2, 'b', tau, x2_0, 'k-.', tau, inter_x2, 'm.');
axis([0, 10, -2, 10]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[-2:2:10]);
grid on;
subplot(3,1,3);
plot(tau, x3, 'b', tau, x3_0, 'k-.', tau, inter_x3, 'm.');
axis([0, 10, 12, 36]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[12:3:36]);
grid on;

figure('Name', 'u');
subplot(2,1,1);
plot(tau, u1, 'b', tau, u1_0, 'k-.', tau, dx1, 'r*');
axis([0, 10, -1, 3]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[-1:0.5:3]);
grid on;
subplot(2,1,2);
plot(tau, u2, 'b', tau, u2_0, 'k-.', tau, dx2, 'r*');
axis([0, 10, -2, 3]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[-2:0.5:3]);
grid on;
