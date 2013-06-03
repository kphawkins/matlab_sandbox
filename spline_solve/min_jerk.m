function [a, ad, add] = min_jerk(x_i, xd_i, xdd_i, x_f, xd_f, xdd_f, d)
A = [  d^3,     d^4,     d^5;
     3*d^2,   4*d^3,   5*d^4;
       6*d,  12*d^2,  20*d^3];

b = [x_f-(x_i+ xd_i*d+xdd_i*d^2);
     xd_f-(xd_i+2*xdd_i*d);
     xdd_f-(2*xdd_i)];

a(1:3) = [x_i, xd_i, xdd_i];
a(4:6) = A\b;
ad = [a(2), 2*a(3), 3*a(4), 4*a(5), 5*a(6)];
add = [2*a(3), 6*a(4), 12*a(5), 20*a(6)];
