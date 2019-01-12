px1 = p(1 : cfg_x.nc);
px2 = p(cfg_x.nc+1 : 2*cfg_x.nc);
px3 = p(2*cfg_x.nc+1 : 3*cfg_x.nc);
pu1 = p(3*cfg_x.nc+1 : 3*cfg_x.nc+cfg_u.nc);
pu2 = p(3*cfg_x.nc+cfg_u.nc+1 : length(p));
x1 = B_x * px1;
x2 = B_x * px2;
x3 = B_x * px3;
u1 = B_u * pu1;
u2 = B_u * pu2;

inter_x1 = zeros(num, 1);
inter_x2 = zeros(num, 1);
inter_x3 = zeros(num, 1);
for i=1:num
    if i == num 
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

figure('Name', 'State');
plot3(x1, x2, x3);
set(gca,'XTick',[10:2:26]);
set(gca,'YTick',[-2:1:8]);
grid on;
figure('Name', 'control input');
plot(u1, u2, dB_x*px1, dB_x*px2, 'r*');
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
plot(tau, u1, 'b', tau, u1_0, 'k-.');
axis([0, 10, -1, 3]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[-1:0.5:3]);
grid on;
subplot(2,1,2);
plot(tau, u2, 'b', tau, u2_0, 'k-.');
axis([0, 10, -2, 3]);
set(gca,'XTick',[0:1:10]);
set(gca,'YTick',[-2:0.5:3]);
grid on;