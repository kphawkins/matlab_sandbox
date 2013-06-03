
clear
% all parameters
%q = [0, 2*pi, pi/2, pi];
%t = [0, 0.5, 2, 3, 4.5, 5];
num_pts = 4;
q = [0, 2*(rand(1,num_pts-2)-0.5)*pi*2, pi];
t = [0, cumsum(2*rand(1,num_pts+1)+0.5)];
qd_i = 0;
qdd_i = 0;
qd_f = 0;
qdd_f = 0;


h = t(2:end)-t(1:end-1);
N = numel(t)-1;
q = [q(1), 0, q(2:end-1), 0, q(end)]; % add in knot placeholders
[a,b,c,d] = make_tri(h,q,qd_i,qdd_i,qd_f,qdd_f); % create tridiagonal matrix problem
X = TDMAsolver(a,b,c,d); % solve tridiagonal matrix problem
QDD = [qdd_i, X, qdd_f]; % create acceleration vector at the knots

for k=1:N
    if k == 1
        C(1) = qd_i + qdd_i*h(1)/2;
    else
        C(k) = QDD(k)*(h(k)+h(k-1))/2 + C(k-1);
    end
    if k == 2
        q(2) = QDD(2)/6*h(1)^2 + C(1)*t(2) + D(1);
    end
    if k == N
        C(k) = qd_f-QDD(k+1)/2*h(k);
        D(k) = q(k+1) - QDD(k+1)/6*h(k)^2 - C(k)*t(k+1);
        q(k) = QDD(k)/6*h(k)^2 + C(k)*t(k+1) + D(k);
    else
        D(k) = q(k) - QDD(k)/6*h(k)^2 - C(k)*t(k);
    end
end

figure(3)
clf
integ = 0;
for k=1:N
    x = linspace(t(k),t(k+1),100);
    y = QDD(k)/(6*h(k))*(t(k+1)-x).^3 + QDD(k+1)/(6*h(k))*(x-t(k)).^3 + C(k)*x + D(k);
    yd = -QDD(k)/(2*h(k))*(t(k+1)-x).^2 + QDD(k+1)/(2*h(k))*(x-t(k)).^2 + C(k);
    ydd = QDD(k)/(h(k))*(t(k+1)-x)+QDD(k+1)/(h(k))*(x-t(k));
    subplot(4,1,4)
    plot(x,y)
    hold on
    plot([t(k),t(k)], [-5,5],'r')
    plot(t, q, 'go')
    subplot(4,1,3)
    plot(x,yd)
    hold on
    plot([t(k),t(k)], [-5,5],'r')
    subplot(4,1,2)
    plot(x,ydd)
    hold on
    plot([t(k),t(k)], [-5,5],'r')
    subplot(4,1,1)
    integ = integ(end) + cumsum(ydd)*(t(k+1)-t(k))/100;
    plot(x,integ)
    hold on
    plot([t(k),t(k)], [-5,5],'r')
end

