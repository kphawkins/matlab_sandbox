function [a, b, c, d] = make_tri(t, q, qd_i, qdd_i, qd_f, qdd_f)
dt_k = t(2) - t(1);
dt_k1 = t(3) - t(2);
a(1) = 0;
b(1) = dt_k/3 + dt_k1/3;
c(1) = dt_k1/6;
d(1) = -q(2)/dt_k + q(1)/dt_k + q(3)/dt_k1 - q(2)/dt_k1 - dt_k/6*qdd_i;
for k = 2:numel(q)-1
dt_k = t(k+1) - t(k);
dt_k1 = t(k+2) - t(k+1);
a(k) = dt_k/6;
b(k) = dt_k/3 + dt_k1/3;
c(k) = dt_k1/6;
d(k) = -q(k+1)/dt_k + q(k)/dt_k + q(k+2)/dt_k1 - q(k+1)/dt_k1;
end
k2 = numel(q);
dt_k = t(k+1) - t(k);
a(k) = dt_k/6;
b(k) = dt_k/3 + dt_k1/3;
c(k) = 0;
d(k) = -q(k+1)/dt_k + q(k)/dt_k + q(k+2)/dt_k1 - q(k+1)/dt_k1 - dt_k1/6*qdd_f;
end
