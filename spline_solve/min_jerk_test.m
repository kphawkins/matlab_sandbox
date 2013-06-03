x_i = rand(1)
x_f = rand(1)
xd_i = 8*(rand(1)-0.5)
xd_f = 8*(rand(1)-0.5)
xdd_i = 0
xdd_f = 0
d = rand(1)+0.5
t = linspace(0,d,1000);
[a,ad,add] = min_jerk(x_i, xd_i, xdd_i, x_f, xd_f, xdd_f, d)
figure(1)
clf
subplot(3,1,3)
plot(t,polyval(a(end:-1:1),t))
hold on
plot([0,d],[x_i,x_f],'go')
subplot(3,1,2)
plot(t,polyval(ad(end:-1:1),t))
hold on
plot([0,d],[xd_i,xd_f],'go')
subplot(3,1,1)
plot(t,polyval(add(end:-1:1),t))
hold on
plot([0,d],[xdd_i,xdd_f],'go')
